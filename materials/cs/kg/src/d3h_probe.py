"""d3h probe — float64/fmpq scale check for the R2-prime Parseval route. NOT certified.

Purpose (minutes): establish magnitudes only.
  (1) sum over grid A^2 == 1 exactly?           (convention lock)
  (2) b_m for odd m in [1,251]: values, P1 = sum m^6 |b_m|^2, head = sum |b_m|
  (3) b_m for odd m in (251,1255]: grid part (exact pairs a+b<=251) — float scale
  (4) nongrid heuristic: (pi/2) sum_a ||g_a||^2 * max_{b} |coeff^a_{m-b}| scale
  (5) bulk mass distribution of m^6 |b_m|^2 vs m — where does the mass live
No claims are certified here.
"""
import json, math, sys, time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gate_A_final import load_grid, A_to_arb  # noqa: E402
import flint  # noqa: E402

KG = Path(__file__).resolve().parent.parent
GRID = KG / "scratch" / "repo" / "grid_M251.json"

# ---- exact rationals for rho coefficients ----
s3 = Fraction(34101124, 10**8)
s5 = Fraction(5276111, 10**8)
V = 1 + s3 * s3 + s5 * s5  # exact
# rho(t) = (t - s3^2 t^3 + s5^2 t^5)/V
# entry tuple (m_int, e10, r_int); A_to_arb returns (mid_arb, half_arb)
import os
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scratch", "repo")
t0 = time.time()
A_dict = load_grid()  # {(a,b): (m_int, e10, r_int)}
entries = [(a, b, t) for (a, b), t in A_dict.items()]
print("n_entries:", len(entries))
meta = None
POW1 = 1 / V          # coeff of t^1 in rho
POW3 = -s3 * s3 / V   # coeff of t^3
POW5 = s5 * s5 / V    # coeff of t^5


def f2(entry) -> tuple[float, float]:
    """entry = (a, b, (m_int, e10, r_int)) -> (mid, radius) of A_{a,b} as floats."""
    a, b, t = entry
    m_int, e10, r_int = t
    scale = 10.0 ** e10
    return int(m_int) * scale, int(r_int) * scale


sq = {}
tot = 0.0
per_slice = {}
for e in entries:
    a, b, t = e
    m, r = f2(e)
    v = m * m  # radius negligible for float64 scale check
    sq[(a, b)] = v
    tot += v
    per_slice[a] = per_slice.get(a, 0.0) + v
print(f"[chk1] sum A^2 = {tot!r}   (want 1; float-sum noise ok)")
print(f"[chk1] max slice mass = {max(per_slice.values())!r}")

# rho^a coefficients cache: dict a -> dict degree k -> Fraction
rho_pow = {0: {0: Fraction(1)}}


def rho_pow_a(a):
    if a in rho_pow:
        return rho_pow[a]
    prev = rho_pow_a(a - 1)
    cur = {}
    for k, c in prev.items():
        for d, cd in ((1, POW1), (3, POW3), (5, POW5)):
            cur[k + d] = cur.get(k + d, Fraction(0)) + c * cd
    rho_pow[a] = cur
    return cur


def coeff(a, k) -> Fraction:
    p = rho_pow_a(a)
    return p.get(k, Fraction(0))


# (2)+(3) b_m over odd m 1..1255 — float scale only. Use exact Fraction numerator
# b_m = (pi/2) sum (−1)^b sq × coeff(a, m−b). Sum in float64 from exact Fractions.
maxdeg_ring = [[] for _ in range(0)]  # unused
bM = {}  # m -> float value (signed)
contrib_count = {}
entries_by_b = {}
for e in entries:
    a, b, _ = e
    entries_by_b.setdefault(b, []).append((a, sq[(a, b)]))
for m in range(1, 1256, 2):
    s = 0.0
    cnt = 0
    for b, lst in entries_by_b.items():
        for a, v in lst:
            k = m - b
            # window a+b<=m<=5a+b ⇔ a<=k<=5a; parity via coeff==0
            if k < a or k > 5 * a:
                continue
            c = coeff(a, k)
            if c == 0:
                continue
            s += ((-1) ** b) * v * float(c)
            cnt += 1
    bM[m] = s
    contrib_count[m] = cnt

b1 = bM[1]
paper_b1 = 0.8815738220495995  # low anchor for scale compare
print(f"[chk3] b1 = {b1!r}  vs paper low {paper_b1!r}  diff {b1 - paper_b1:.3e}")
bad_even = []
for m in range(2, 1256, 2):
    for b, lst in entries_by_b.items():
        for a, v in lst:
            k = m - b
            if k >= a and k <= 5 * a and (k - a) % 2 == 0:
                bad_even.append((m, a, b, k))
print(f"[chk2] even-m nonzero-count (should be 0): {len(bad_even)}")

# head sum & P1 & mass distribution
head = sum(abs(bM[m]) for m in range(3, 252, 2))
P1 = sum((m**6) * bM[m] ** 2 for m in range(1, 252, 2))
print(f"[chk4] sum_3..251 |b_m|      = {head:.6e}   (paper: 1.13288599276897e-5)")
print(f"[chk4] P1 (m<=251) m^6|b|^2  = {P1:.6e}   sqrt = {math.sqrt(P1):.6f}")

grid_tail = 0.0
Pg = 0.0
for m in range(253, 1256, 2):
    grid_tail += abs(bM[m])
    Pg += m**6 * bM[m] ** 2
print(f"[chk5] grid-part (251,1255] sum|b_m| = {grid_tail:.3e};  Σ m^6 b^2 = {Pg:.3e}")

# (4) nongrid heuristic: need ||g_a||^2 = E[q_a(ϑψ3(X))^2] over those a with
# a + b > 251 achievable for some odd m ≤ 1255 band; quick scale: use ||g_a||^2 ~ c_a from feed
# but cheap estimate: for each a present use per_slice mass of a IF a<=251-like else ≤1 heuristic.
vartheta = float(Fraction(136419125, 10**9)) / math.sqrt(float(V))


def He(n, x):
    """Probabilists' Hermite He_n via forward recurrence He_{n+1}=x He_n − n He_{n−1}."""
    h_prev, h_cur = 1.0, x
    if n == 0:
        return h_prev
    if n == 1:
        return h_cur
    for k in range(1, n):
        h_prev, h_cur = h_cur, x * h_cur - k * h_prev
    return h_cur


def q_a(a, u):
    """Vectorized over u (numpy). q_0 = erf(u/√2); q_a = 2φ(u)He_{a−1}(−u)/√a."""
    u = np.asarray(u)
    if a == 0:
        return np.vectorize(math.erf)(u / math.sqrt(2.0))
    # He_{a-1}(−u) by recurrence on the array
    n = a - 1
    hp = np.ones_like(u)
    if n == 0:
        hc = hp
    else:
        hc = -u  # He_1(x) = x evaluated at x=−u
        for k in range(1, n):
            hp, hc = hc, (-u) * hc - k * hp
    phi_u = np.exp(-u * u / 2.0) / math.sqrt(2 * math.pi)
    return 2.0 * phi_u * hc / math.sqrt(a)


# slice masses ||g_a||^2 = E[q_a(ϑ ψ3(X))^2] via coarse Gauss–Hermite float
import numpy as np  # noqa: E402

x, w = np.polynomial.hermite_e.hermegauss(220)  # nodes for weight e^{−x²/2}
sl = {}
for a in range(0, 400):
    xc = np.clip(x, -7.0, 7.0)  # |u(7)|≈17; He_{a−1}(17) ≤ 17^99 ≈ 1e122 OK; farther nodes carry < 1e-40 mass
    u = vartheta * (xc**3 - 3.0 * xc) / math.sqrt(6.0)
    val = q_a(a, u)
    sl[a] = float(np.dot(np.clip(w, 0, None), val * val) / math.sqrt(2.0 * math.pi))
print(f"[chk6] slice masses: s_1={sl[1]:.4f} s_3={sl[3]:.4f} s_10={sl[10]:.5f} s_50={sl[50]:.3e} "
      f"s_{max(sl)}={sl[max(sl)]:.3e}, total(0..399)={sum(sl.values()):.6f}")

# nongrid mass heuristic on (251,1255]: for each m, Σ_a [mass beyond b>251−a in slice a × coeff bound],
# using slice mass ≤ min(1, sl[a] if a<400 else est) — coarse.
# tau_a(m): max coeff |[t^k]ρ^a| over achievable k=a..5a for given m-b range...
# Cheap screening first: total possible nongrid mass Σ over a of (1-per_slice-mass)*(bandwidth).
# Actual estimate: for each odd m in band, nongrid ≤ (π/2)*Σ_a sl.get(a,≤1)*M_a(m),
# M_a(m) = max_k in [a,5a], k≡a (2), b=m−k>251−a → k ranges where m−k>251−a ⇔ k<m−251+a.
Mn = 0.0   # max of (pi/2)*Σ_a clause over m
tail_m_estimators = {}
for m in range(253, 1256, 2):
    sN = 0.0
    for a in range(1, m):
        # nongrid pairs (a,b): b ≥ 252−a (integer b > 251−a), k = m−b, need a ≤ k ≤ 5a, k ≡ a (mod 2)
        hi = min(5 * a, m - 252 + a)
        lo = a
        if lo > hi:
            continue
        lo2 = lo if (lo - a) % 2 == 0 else lo + 1
        if lo2 > hi:
            continue
        # probe k values maximizing |coeff| among parity-matched (coarse)
        cands = [lo2, hi - ((hi - a) % 2), (lo2 + hi) // 2]
        best = 0.0
        for k in cands:
            if a <= k <= 5 * a and (k - a) % 2 == 0:
                best = max(best, abs(float(coeff(a, k))))
        sN += min(sl.get(a, 1.0), 1.0) * best
    val = (math.pi / 2) * sN
    tail_m_estimators[m] = val
    Mn = max(Mn, val)
est_sum = sum(tail_m_estimators.values())
est_P2 = sum((m**6) * tail_m_estimators[m] ** 2 for m in range(253, 1256, 2))
print(f"[chk7] nongrid heuristic: max_m=(π/2)Σ_a.. = {Mn:.3e}; Σ|b_nongrid| est = {est_sum:.3e}; "
      f"Σ m^6 b^2 (nongrid, OVEREST since squared) = {est_P2:.3e}")

# combined sqrt estimate
totP = P1 + Pg + est_P2
print(f"[chk8] sqrt(P1+Pg+estP2) = {math.sqrt(totP):.6f}   (authors: 14.44243664663976457; (π/√6)√126.804 = {math.sqrt(math.pi**2/6*126.80385221):.4f})")
print(f"[probe done in {time.time()-t0:.1f}s]")
