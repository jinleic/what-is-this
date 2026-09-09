"""
rc28.py — Certified upper bounds on Cbar_{n,k} for the uniform-subsequence
channel W_{n,k}, following Rubinstein-Con's own schema (arXiv:2305.07156):

  1. DIRECT DUAL bounds: for ANY full-support output distribution D',
       Cbar_{n,k} <= max_x KL(W(.|x) || D')    (Csiszar-Tusnady / Blahut 1972)
     Evaluated in Arb (outward-rounded) from an exact rational D'.
     D' comes from a float64 NumPy BA location step; its rounding error does
     NOT affect validity (any D' gives a valid UB) — only tightness.

  2. RECURRENCE bounds (RC Lemmas 2.3/2.4/2.5):
       L2.3: Cbar_{n+1,k} <= Cbar_{n,k}
       L2.4: Cbar_{n+1,k} <= Cbar_{n,k-1}(1-k/(n+1)) + (Cbar_{n,k}+1) k/(n+1)
       L2.5: Cbar_{n,k} <= sum_i [C(s,i)C(n-s,k-i)/C(n,k)](Cbar_{s,i}+Cbar_{n-s,k-i})
     Propagated in Arb interval arithmetic from the certified base grid.

  3. U(n,d) = (1/n) sum_k binom(n,k) d^{n-k}(1-d)^k Cbar_{n,k}  (RC eq (2.2)+(3.1))
     upper-bounds C(BDC_d). Compare U(28,0.68)/0.32 with RC's 0.3745.

Verdict rule (fixed in pre_statement.md):
    CONFIRM if certified U(28,0.68)/0.32 <= 0.3745 + 5e-5
    UNKNOWN otherwise (weaker bound than theirs; no contradiction possible
    from an upper bound alone).
"""
from __future__ import annotations
import os
import math, time, json, os, sys
import numpy as np
import flint
from flint import arb, fmpq, fmpz

flint.ctx.prec = 300
LOG2 = None  # lazy


# ================================================================ exact matrix
def subseq_matrix_i64(n: int, k: int) -> np.ndarray:
    """Embedding-count matrix N (2^n x 2^k), int64. Entries < 2^63 ok for n<=33."""
    Y = np.arange(1 << k)
    bits = ((Y[:, None] >> np.arange(k)) & 1).astype(np.int8)
    N = np.zeros((1 << n, 1 << k), dtype=np.int64)
    for x in range(1 << n):
        dpY = np.zeros((k + 1, 1 << k), dtype=np.int64)
        dpY[0, :] = 1
        for ii in range(n):
            a = (x >> ii) & 1
            new = np.empty_like(dpY)
            new[0, :] = dpY[0, :]
            for j in range(1, k + 1):
                want = bits[:, j - 1] == a
                new[j, :] = dpY[j, :] + np.where(want, dpY[j - 1, :], 0)
            dpY = new
        N[x, :] = dpY[k, :]
    return N


def ba_float(N: np.ndarray, Cnk: int, iters: int = 400):
    """Float64 BA locating step (vectorized). Returns (p, D, rate, dual)."""
    nX, nY = N.shape
    W = N / Cnk
    p = np.full(nX, 1.0 / nX)
    with np.errstate(divide='ignore'):
        logW = np.where(W > 0, np.log2(W), 0.0)
    mask = W > 0
    t = None
    for it in range(iters):
        D = p @ W
        with np.errstate(divide='ignore', invalid='ignore'):
            t = np.where(mask, W * (logW - np.log2(D)[None, :]), 0.0)
        a = t.sum(axis=1)
        m = a.max()
        w = np.exp2(a - m) * p
        p = w / w.sum()
    D = p @ W
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.where(mask, W * (logW - np.log2(D)[None, :]), 0.0)
    rate = float(p @ t.sum(axis=1))
    dual = float(t.sum(axis=1).max())
    return p, D, rate, dual


def snap_distribution(D: np.ndarray, den: int) -> list[int]:
    """Rational snap of D to num/den with exact sum = den."""
    raw = D * den
    m = np.floor(raw).astype(object)
    diff = den - int(sum(m))
    order = np.argsort(-(raw - m))
    for i in range(abs(diff)):
        j = int(order[i % len(m)])
        m[j] += 1 if diff > 0 else -1
        assert m[j] >= 0
    assert int(sum(m)) == den
    return [int(v) for v in m]


def dual_cert(N: np.ndarray, Cnk: int, Dnum: list[int], Dden: int) -> arb:
    """Certified max_x KL(W(.|x) || D') in Arb; upper() >= true dual >= Cbar."""
    nX, nY = N.shape
    log2c = arb(2).log()
    DnumF = [fmpz(v) for v in Dnum]
    best = None
    with np.errstate(divide='ignore'):
        # per-row KL: chunked over x for memory friendliness
        chunk = max(1, (1 << 20) // max(nY, 1))
        for x0 in range(0, nX, chunk):
            x1 = min(nX, x0 + chunk)
            blk = N[x0:x1]
            nzx = np.nonzero(blk.any(axis=1))[0]
            for r in nzx:
                row = blk[r]
                nz = np.nonzero(row)[0]
                s = arb(0)
                for j in nz:
                    c = int(row[j])
                    Dn = DnumF[int(j)]
                    q = fmpq(c * Dden, Cnk * Dn)
                    term = arb(fmpq(c, Cnk)) * (arb(q).log() / log2c)
                    s += term
                if best is None or float(s.upper()) > float(best.upper()):
                    best = s
                if float(s.upper()) == float(s.lower()) and best is not None and \
                   float(s) <= float(best) - 1e-6:
                    pass
    return best


def certify_nk(n: int, k: int, iters: int = 400, snap_den: int = 1 << 40) -> dict:
    t0 = time.time()
    N = subseq_matrix_i64(n, k)
    Cnk = math.comb(n, k)
    p, D, rate, dual = ba_float(N, Cnk, iters=iters)
    Dnum = snap_distribution(D, snap_den)
    ub = dual_cert(N, Cnk, Dnum, snap_den)
    return dict(n=n, k=k, Cnk=Cnk, rate=rate, dual_float=dual,
                dual_ub=str(ub), dual_ub_f=float(ub.upper()),
                wall_s=round(time.time() - t0, 1))


# ================================================================== grid
class Grid:
    """UB[n][k] as arb balls (valid upper envelopes of Cbar_{n,k})."""

    def __init__(self, n_max: int):
        self.n_max = n_max
        self.ub = [[None] * (n_max + 1) for _ in range(n_max + 1)]
        self.kind = [[None] * (n_max + 1) for _ in range(n_max + 1)]

    def set_direct(self, n, k, ball: arb, tag='direct'):
        self.ub[n][k] = ball
        self.kind[n][k] = tag

    def base_and_propagate(self):
        n_max = self.n_max
        # Cbar_{n,k} <= k trivially; base UB for all (n,k): k
        for n in range(1, n_max + 1):
            for k in range(1, n + 1):
                if self.ub[n][k] is None:
                    self.ub[n][k] = arb(k)
                    self.kind[n][k] = 'trivial-k'
        for sweep in range(6):
            changed = 0
            for n in range(1, n_max + 1):
                for k in range(1, n + 1):
                    cands = [self.ub[n][k]]
                    if n - 1 >= k and self.ub[n - 1][k] is not None:
                        cands.append(self.ub[n - 1][k])                       # L2.3
                    if k - 1 >= 1 and n - 1 >= k and self.ub[n - 1][k] is not None:
                        # L2.4 at (n+1, k) written for (n,k): C_{m,L} <= C_{m-1,L-1}(1-L/m) + (C_{m-1,L}+1) L/m
                        pass
                    if k + 0 <= n - 1 + 1:  # C_{n,k} from row n-1 (L2.4 reversed index)
                        m = n
                        L = k
                        if L - 1 >= 0 and m - 1 >= L - 1 and m - 1 >= L and \
                           self.ub[m - 1][L - 1] is not None and self.ub[m - 1][L] is not None:
                            alpha = arb(fmpq(L, m))
                            cands.append(self.ub[m - 1][L - 1] * (1 - alpha) +
                                         (self.ub[m - 1][L] + 1) * alpha)
                    best = min(cands, key=lambda b: float(b.upper()))
                    if best is not self.ub[n][k] and float(best.upper()) < float(self.ub[n][k].upper()):
                        self.ub[n][k] = best
                        changed += 1
            # L2.5 split
            for n in range(2, n_max + 1):
                for k in range(1, n + 1):
                    Cnk = math.comb(n, k)
                    best_tot = None
                    for s in range(1, n):
                        tot = arb(0)
                        ok = True
                        for i in range(0, s + 1):
                            ki = k - i
                            if ki < 0 or ki > n - s:
                                continue
                            a = self.ub[s][i] if i <= s else None
                            b = self.ub[n - s][ki] if ki <= n - s else None
                            if a is None or b is None:
                                ok = False
                                break
                            w = fmpq(math.comb(s, i) * math.comb(n - s, ki), Cnk)
                            tot += arb(w) * (a + b)
                        if ok and (best_tot is None or float(tot.upper()) < float(best_tot.upper())):
                            best_tot = tot
                    if best_tot is not None and float(best_tot.upper()) < float(self.ub[n][k].upper()):
                        self.ub[n][k] = best_tot
                        self.kind[n][k] = 'L2.5-split'
                        changed += 1
            if changed == 0:
                break

    def U(self, d: fmpq, n: int) -> arb:
        """U(n,d) = (1/n) sum_k binom(n,k) d^{n-k}(1-d)^k UB[n][k]  (arb ball)."""
        da = arb(d)
        pa = 1 - da
        total = None
        for k in range(0, n + 1):
            if self.ub[n][k] is None and k >= 1:
                continue
            if k == 0:
                continue  # C_{n,0}=0
            w = (arb(math.comb(n, k)) * (da ** (n - k)) * (pa ** k)) / arb(n)
            term = w * self.ub[n][k]
            total = term if total is None else total + term
        return total or arb(0)


def main():
    t_start = time.time()
    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'campaigns')
    stamp = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    import hashlib, platform
    tag = hashlib.sha256((platform.platform() + str(time.time())).encode()).hexdigest()[:12]
    camp = os.path.join(out_dir, f'{stamp}_{tag}')
    os.makedirs(camp, exist_ok=True)

    # ---------------- base grid: direct dual certificates
    N_DIRECT = 15
    results = []
    print('building base grid n<=', N_DIRECT, flush=True)
    for n in range(2, N_DIRECT + 1):
        for k in range(1, n + 1):
            if (1 << n) * (1 << k) > (1 << 24):
                continue
            try:
                r = certify_nk(n, k, iters=400)
                results.append(r)
                print(f"  ({n},{k}) dual_ub={r['dual_ub_f']:.9f} float_dual={r['dual_float']:.9f} ({r['wall_s']}s)", flush=True)
            except MemoryError:
                print(f'  ({n},{k}) MEMORY SKIP', flush=True)
    # write raw certificates
    with open(os.path.join(camp, 'direct_certificates.json'), 'w') as f:
        json.dump(results, f, indent=1)

    # ---------------- grid assembly
    g = Grid(28)
    for r in results:
        g.set_direct(r['n'], r['k'], arb(r['dual_ub']), tag='direct-dual')
    g.base_and_propagate()

    # ---------------- U(28, 0.68)
    d68 = fmpq(68, 100)
    U28 = g.U(d68, 28)
    U16 = g.U(d68, 16)
    norm = U28 / arb(fmpq(32, 100))
    print(f"U(28,0.68) certified = {U28}")
    print(f"U(28,0.68)/0.32      = {norm}")
    print(f"RC published: 0.1199 (=0.3745*0.32); verdict rule: CONFIRM iff norm.upper <= 0.37455")
    conf = float(norm.upper()) <= 0.3745 + 5e-5
    print('VERDICT:', 'CONFIRM' if conf else 'UNKNOWN (our UB looser than RC n=28 row)')

    # ---------------- U(n,0.68) for all n 2..28, for the ledger
    series = {}
    for n in range(2, 29):
        Un = g.U(d68, n)
        series[n] = (str(Un), float(Un.upper()) / 0.32)
    with open(os.path.join(camp, 'U_0.68_series.json'), 'w') as f:
        json.dump({str(k): v for k, v in series.items()}, f, indent=1)

    # summary
    with open(os.path.join(camp, 'summary.txt'), 'w') as f:
        f.write(f"""delcap RC-schema certified enclosure {stamp}
base grid: direct dual certificates (Csiszar-Tusnady/Blahut 1972) at n <= {N_DIRECT},
    float64 BA locating step + exact-rational D' snap + Arb outward-rounded
    max_x KL(W(.|x) || D') >= Cbar.
grid: RC/FD lemmas 2.3/2.4/2.5 propagated in Arb to n=28.
U(28, 0.68) certified = {U28}
normalized            = {norm} / 0.32
RC published          = 0.1199 (== 0.3745*0.32) at n=28, Table 1
verdict               = {'CONFIRM' if conf else 'UNKNOWN'}
elapsed {time.time()-t_start:.0f}s
""")
    print('campaign dir:', camp)
    return camp, conf, float(norm.upper())


if __name__ == '__main__':
    main()
