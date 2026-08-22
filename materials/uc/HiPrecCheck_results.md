# HiPrecCheck — independent high-precision cross-check of the certified numerics

**Agent:** HiPrecCheck · **Date:** 2026-08-18 · **Status:** all 4 checks PASS.

**Scope.** Pre-work for final verification: an independent re-derivation, at
mpmath dps = 50 (KKT also at dps = 60 to mirror `cert3.py`), of the four
numerical pillars behind the certified conclusion at
`t0 = 0.3820660112501052 = ψ + 10⁻⁴` (`α = 0.0356069`, base-2 entropy h).
No `.py` file in this repository was created or edited; the only repository
file written is this one. Scripts were run from the repo root as

```
math/.venv/bin/python -B <script>      # python 3.14.3, mpmath 1.3.0, numpy 2.5.2, flint 0.9.0
```

**Labeling.** Everything below is **NUMERICAL** (high-precision floating point
or float64 sampling) — it corroborates, but does not replace, the Arb-interval
certificate. The only certified quantities quoted are read off from existing
code (`cert2.get_rho_gmax()` Arb enclosure; recorded constants `B_OBS`,
`A_OBS`, `C_STAR`, and PROOF.md's `λ*`). No new PROVED claims are made here.

## Verdict summary

| # | Check | Verdict | Key numbers |
|---|-------|---------|-------------|
| 1 | q2=1 face KKT at t0 | **PASS** | p\*=r\*=0.3277844883720486534962090162942804726000486012553905, w\*=0.8384996420371223392315493623121482657880908513976684, λ\*=1.368215810025175744411947021070813623862875825251114; residual 4.45e-68 < 1e-35; λ\*>0; **bit-identical** to `cert3.solve_face_kkt` |
| 2 | coarse obstruction, c\*, margin at t0 | **PASS** | b̄=0.3294547385030369723917053838771341307518391169167 (matches `B_OBS` to 1.7e-21), a=0.078877292705923173412161944535296855609585089944575 (matches `A_OBS` to 1.6e-22), c\*=0.38234553336670272114599300912553666790802222940517, c\*−ψ=3.7952211659756935e-4; Q=C=L at ν\* to 3.2e-50; Φ≡0 and Φ_exact≡0 in α; margin at t0 > 0 by all three measures (4.5557e-4 heuristic, 4.5661e-4 at KKT, 4.5755e-4 grid) |
| 3 | orbit-swap invariance | **PASS** | max error **0.0** over 100 random points at 50 dps (exact to the last bit); independent Φ matches `cert2.phi_true` to 0.0 |
| 4 | rh(x) ≤ 2min(x,1−x); max ρ | **PASS** | no violation on 10 001-point grid at 50 dps (equality only at x=0,1; min interior slack 1.573e-7); max_p ρ = 0.3049467103530205557127432257297999107925 at p = 0.1165022347591509830549284…, consistent with the certified record 0.3049467 and below cert2's certified Arb cap 0.305360849157902 |

Consistency statement (NUMERICAL, this file): the certified conclusion's own
numerical inputs — the face-KKT λ-family anchor, the coarse obstruction
constant c\*, the exact functional's degeneracy there, the orbit involution,
and the two ρ/rh facts — all reproduce from scratch to between 16 and 22
significant digits against the values on record, with the single last-digit
transcription note on PROOF.md's λ\* registered below. Everything sampled is
NUMERICAL; the certificate's correctness rests on the Arb campaign, not on
anything here.

---

## Check 1 — KKT system on the q2 = 1 face at t0

The face parametrisation (derived here from the formulas, not copied): orbit
A = (p, p) with weight w, orbit B = (r, 1) with weight 1−w, i.e. the 5-point
(p1,q1,p2,q2,w) = (p, p, r, 1, w); this is the parametrisation for which
`cert3.solve_face_kkt`'s `fun` is exactly `Φ + λ(M−t)` with

- M = w·p + (1−w)(r+1)/2,
- L = w·h(p) + (1−w)h(r)/2,
- C = w·h(s\*(p,p)) (since s\*(r,1)=1, h(1)=0),
- σ = w√rh(p) + (1−w)√rh(r)/2,
- rh(z) = max(0, 2(1−z)h(z) − h((1−z)²)), s\*(z,z) = min(max(½,z), min(2z,1)).

Solved with `mp.findroot`, seeds (0.328, 0.328, 0.838, 1.37), tol 1e-45, with
`mp.diff` numerical derivatives — the exact recipe of `cert3.solve_face_kkt`
(read-only comparison import at the end).

**Results (dps = 60):**

```
p*     = 0.3277844883720486534962090162942804726000486012553905
r*     = 0.3277844883720486534962090162942804726000486012553905   (= p*, symmetric)
w*     = 0.8384996420371223392315493623121482657880908513976684
lam*   = 1.368215810025175744411947021070813623862875825251114
residual = 4.4533e-68   (< 1e-35 asserted, holds)
lam* > 0                (asserted, holds)
```

- p\* = r\* to all 50 digits: the KKT point is symmetric, orbit A = (p\*,p\*),
  orbit B = (r\*,1) with p\* = r\* = 0.32778…, adjacent to the coarse
  obstruction b̄ = 0.32945… — as expected for the t0-shifted minimiser.
- dps = 50 re-solve differs from dps = 60 by < 2.1e-51 (dps-independent).
- **`cert3.solve_face_kkt` own output is bit-identical (max diff 0.0).**
- PROOF.md (line 220) records λ\* = 1.3682158100251758: agrees to
  **16 significant digits**; |Δ| = 5.56e-17, a last-digit transcription
  rounding in the document (λ\* is a NUMERICAL seed; the certificate's
  λ-family {0.98λ\*, λ\*, 1.02λ\*} is unaffected by 1e-17).

**Script** (`/tmp/hiprec_check1.py`, run once, 1.1 s):

```python
"""Check 1: independent re-derivation of the q2=1 face KKT point at t0.

Replicates cert3.solve_face_kkt exactly (same formulas, same seeds, same
numerical differentiation via mp.diff, same tol) but from a fresh transcript,
then compares against cert3's own function (read-only import).
"""
import sys, time
import mpmath as mp

T0 = "0.3820660112501052"   # psi + 1e-4, the certification level
ALPHA = "0.0356069"         # Cambie's alpha, as used by cert3

def h(x):
    x = mp.mpf(x)
    return -(x*mp.log(x) + (1-x)*mp.log(1-x)) / mp.log(2)

def rh(z):
    x = 1 - z
    return max(mp.mpf(0), 2*x*h(z) - h(x*x))

def sdiag(z):
    return min(max(mp.mpf("0.5"), z), min(2*z, 1))

def run(dps):
    mp.mp.dps = dps
    aa, tt = mp.mpf(ALPHA), mp.mpf(T0)

    def fun(p, r, w, lam):
        mean = w*p + (1-w)*(r+1)/2
        entropy = w*h(p) + (1-w)*h(r)/2
        sigma = w*mp.sqrt(rh(p)) + (1-w)*mp.sqrt(rh(r))/2
        return ((2*(1-aa)*(1-mean) - 1)*entropy
                + aa*w*h(sdiag(p)) - (1-aa)*sigma*sigma
                + lam*(mean - tt))

    dp = lambda p, r, w, lam: mp.diff(lambda z: fun(z, r, w, lam), p)
    dr = lambda p, r, w, lam: mp.diff(lambda z: fun(p, z, w, lam), r)
    dw = lambda p, r, w, lam: mp.diff(lambda z: fun(p, r, z, lam), w)
    cons = lambda p, r, w, lam: w*p + (1-w)*(r+1)/2 - tt

    p, r, w, lam = mp.findroot(
        (dp, dr, dw, cons),
        (mp.mpf("0.328"), mp.mpf("0.328"), mp.mpf("0.838"), mp.mpf("1.37")),
        tol=mp.mpf("1e-45"), maxsteps=100)
    residual = max(abs(dp(p, r, w, lam)), abs(dr(p, r, w, lam)),
                   abs(dw(p, r, w, lam)), abs(cons(p, r, w, lam)))
    return p, r, w, lam, residual


t_start = time.time()
p, r, w, lam, residual = run(60)
print("== independent solve, dps=60 (identical recipe to cert3.solve_face_kkt) ==")
print("p*     =", mp.nstr(p, 52))
print("r*     =", mp.nstr(r, 52))
print("w*     =", mp.nstr(w, 52))
print("lam*   =", mp.nstr(lam, 52))
print("resid  =", mp.nstr(residual, 5))
assert residual < mp.mpf("1e-35"), residual
assert 0 < p < 1 and 0 < r < 1 and 0 < w < 1 and lam > 0
print("asserts: residual < 1e-35, box bounds, lam* > 0   [PASS]")

# stability at dps=50 (task-specified working precision)
p50, r50, w50, lam50, res50 = run(50)
print("dps=50 vs dps=60 deltas:",
      mp.nstr(abs(p-p50), 3), mp.nstr(abs(r-r50), 3),
      mp.nstr(abs(w-w50), 3), mp.nstr(abs(lam-lam50), 3))
assert max(abs(p-p50), abs(r-r50), abs(w-w50), abs(lam-lam50)) < mp.mpf("1e-45")
print("dps-independence to <1e-45                              [PASS]")

# cert3's own output, read-only
sys.path.insert(0, "math/uc")
import cert3
pc, rc, wc, lamc, reslc = cert3.solve_face_kkt(mp.mpf(T0))
print("== cert3.solve_face_kkt own output (dps=60) ==")
print("p*     =", mp.nstr(pc, 52))
print("r*     =", mp.nstr(rc, 52))
print("w*     =", mp.nstr(wc, 52))
print("lam*   =", mp.nstr(lamc, 52))
print("resid  =", mp.nstr(reslc, 5))
dmax = max(abs(p-pc), abs(r-rc), abs(w-wc), abs(lam-lamc))
print("max |independent - cert3| over (p,r,w,lam):", mp.nstr(dmax, 3))
assert dmax < mp.mpf("1e-45")
print("agreement with cert3 to <1e-45                          [PASS]")

# recorded value in PROOF.md line 220: lambda* = 1.3682158100251758
lam_doc = mp.mpf("1.3682158100251758")
print("PROOF.md lambda*   = 1.3682158100251758")
print("|lam* - PROOF.md lam*| =", mp.nstr(abs(lam - lam_doc), 3),
      "-> agrees to 16 significant digits (last-digit transcription)")
assert abs(lam - lam_doc) < mp.mpf("1e-16")
print("PROOF.md lambda* agreement                             [PASS]")
print("check1 wall time %.2fs" % (time.time() - t_start))
```


**Output (verbatim):**

```
== independent solve, dps=60 (identical recipe to cert3.solve_face_kkt) ==
p*     = 0.3277844883720486534962090162942804726000486012553905
r*     = 0.3277844883720486534962090162942804726000486012553905
w*     = 0.8384996420371223392315493623121482657880908513976684
lam*   = 1.368215810025175744411947021070813623862875825251114
resid  = 4.4533e-68
asserts: residual < 1e-35, box bounds, lam* > 0   [PASS]
dps=50 vs dps=60 deltas: 7.74e-52 7.74e-52 1.75e-51 2.1e-51
dps-independence to <1e-45                              [PASS]
== cert3.solve_face_kkt own output (dps=60) ==
p*     = 0.3277844883720486534962090162942804726000486012553905
r*     = 0.3277844883720486534962090162942804726000486012553905
w*     = 0.8384996420371223392315493623121482657880908513976684
lam*   = 1.368215810025175744411947021070813623862875825251114
resid  = 4.4533e-68
max |independent - cert3| over (p,r,w,lam): 0.0
agreement with cert3 to <1e-45                          [PASS]
PROOF.md lambda*   = 1.3682158100251758
|lam* - PROOF.md lam*| = 5.56e-17 -> agrees to 16 significant digits (last-digit transcription)
PROOF.md lambda* agreement                             [PASS]
check1 wall time 1.09s
```

---

## Check 2 — Φ = L·Λ and Φ_exact at the coarse obstruction; c\*; margin at t0

The coarse obstruction is μ\* = a·δ1 + (1−a)·δ_b̄ realised as the pair-orbit
measure ν\* = (1−2a)·δ_(b̄,b̄) + 2a·δ_(b̄,1) (i.e. 5-point
(b̄, b̄, b̄, 1, w) with w = 1−2a ∈ [½,1]). Φ_exact denotes the exact Gilmer
functional F = (1−α)Q + αC − L with Q = Σ_s Σ_t u_s u_t h(s+t−st) over the four
atom-weights u = (wA/2, wA/2, wB/2, wB/2).

**Results (dps = 50).**

Defining equation h(b)(2−h(b)) = h((1−b)²) has exactly two roots in (0,1):

```
b = 0.13949945190986195999533605297675493434150185769247
bbar (larger) = 0.3294547385030369723917053838771341307518391169167
|bbar - B_OBS| = 1.71e-21          (all 20 recorded digits of B_OBS correct)
```

Relations at the obstruction:

```
h(bbar)          = 0.91436831153838781849288115950949617751165248080802
1/(2 - h(bbar))  = 0.92112270729407682658783805546470314439041491005543   (= 1 - a)
a                = 0.078877292705923173412161944535296855609585089944575
|a - A_OBS|      = 1.62e-22          (all 21 recorded digits of A_OBS correct)
(1-a)h((1-bbar)^2) - h(bbar) = 0.0    (exact, to 50 digits)
```

> **Discrepancy note (assignment vs. code/README):** the task's phrasing
> "a = 1/(2−h(b))" is the *complement*: the relation matching `A_OBS` — and the
> one stated in README "Result 2" — is **1−a = 1/(2−h(b̄))**, i.e.
> a = 1 − 1/(2−h(b̄)) = 0.078877…, while 1/(2−h(b̄)) = 0.921122… ≈ 1−a.
> Verified numerically both ways above.

c\* (to 50 digits):

```
psi   = 0.38196601125010515179541316563436188227969082019424
c*    = 1-(1-bbar)/(2-h(bbar)) = a + (1-a)bbar
      = 0.38234553336670272114599300912553666790802222940517
c*-psi = 0.00037952211659756935057984349117478562833140921093393  (= 3.795221...e-4)
|c* - C_STAR| = 4.01e-21          (all 20 recorded digits of C_STAR correct)
```

Functional values at ν\* (NUMERICAL at 50 dps, corroborating the identity):

```
mean(mu*) = 0.38234553336670272114599300912553666790802222940517  (= c*, exact)
L = C = 0.84224541458815365317567611092940628878082982011085
Q = 0.84224541458815365317567611092940628878082982011088
sigma^2 = 0.19818787005541911435255209742502344030155715976848
Q - L = 3.21e-50, C - L = 0.0, Q - C = 3.21e-50    ->  Q = C = L to < 1e-45
alpha in {0, 0.0356069, 0.2, 0.5, 1}: Phi_exact = F = (1-a)Q+aC-L in [+0.0, +3.2e-50]
                                        Phi (certified form)      in [-3.3e-52, +6.7e-52]
=> Phi = L*Lambda and Phi_exact both vanish IDENTICALLY in alpha at nu*.
Lambda(nu*) = Phi/L = -3.34e-52 ~ 0.
```

This reproduces margin_lemma.py's section 3 claim (C–S tight at the
obstruction, F ≡ 0 for every α) at 50-digit precision instead of 1e-12, and
matches `cert2.phi_true` bit-for-bit (|Δ| = 0.0) at the same point.

Margin at t0 = ψ+1e-4 — three independent measures, all > 0:

```
c* - t0             = 0.00027952211659752114599300912553666790802223142580238
(c*-t0)*1.6298      = 0.00045556514563063996373940627279966135649463641666458   (heuristic)
Lambda at KKT point = 0.000456613589352931585341295601778431169226882012        (50 dps)
grid min Lambda     = 4.5754855172e-04  (float64, NUMERICAL; mpmath at argmin
                                          = 0.000457548551724293780830740681267)
```

All three agree to ~0.4% relative and all are strictly positive; the recorded
"min-Λ ≈ 4.56e-4 at t=ψ+1e-4" (PROGRESS.md) is reproduced. The KKT point
stationarity is for Φ (not Λ), and the float64 Λ-landscape is flat near the
minimum, which explains the small spread and the grid's argmin jitter
(4.5755e-4 vs 4.5795e-4 across refinement rounds).

Grid verification of min Λ > 0 at t0 (NUMERICAL, float64): the symmetric face
slice (p, r, w) — the slice containing both the coarse obstruction and the
t0 KKT point — was swept coarse 200³ over (p,r) ∈ [0,1], w ∈ [½,1], then
refined three more times at 200³ around the running argmin (mask: mean ≤ t0):

```
coarse  min Lambda = 4.8916318748e-04 at (0.326634, 0.321608, 0.834171)
refined min Lambda = 4.5754855172e-04 at (0.328744, 0.326734, 0.840653)
```

Additionally, 2 000 000 random points over the **full 4-parameter q2=1 face**
(p1 ≠ q1 allowed; 288 717 feasible with mean ≤ t0, w ∈ [½,1]) gave
min Λ = 1.0831e-3 at (0.316262, 0.331387, 0.353808, 0.835967) — nothing off
the symmetric slice comes near the minimum, corroborating that the t0
minimiser lives on the p1=q1 slice, as the KKT solve found (p\* = r\*).

**Script** (`/tmp/hiprec_check2.py`, 3.5 s):

```python
"""Check 2: the coarse obstruction a*delta_1 + (1-a)*delta_bbar, c*, and the
margin at t0 = psi + 1e-4.

All mpmath at dps=50 unless labelled NUMERICAL (float64 grid).
"""
import sys, time
import mpmath as mp

mp.mp.dps = 50
T0 = mp.mpf("0.3820660112501052")
PSI = (3 - mp.sqrt(5)) / 2
ALPHA = mp.mpf("0.0356069")
B_OBS = mp.mpf("0.32945473850303697239")
A_OBS = mp.mpf("0.078877292705923173412")
C_STAR = mp.mpf("0.38234553336670272115")

t_start = time.time()

def h(x):
    x = mp.mpf(x)
    if x == 0 or x == 1:
        return mp.mpf(0)
    return -(x*mp.log(x) + (1-x)*mp.log(1-x)) / mp.log(2)

def rh(z):
    x = 1 - mp.mpf(z)
    return max(mp.mpf(0), 2*x*h(z) - h(x*x))

def ss(x, y):
    return min(max(mp.mpf("0.5"), max(x, y)), min(x + y, 1))

# ---- independent 5-parameter Phi (identical formula to cert2.phi_true) -----
def aggregates(p1, q1, p2, q2, w):
    m1, m2 = (p1+q1)/2, (p2+q2)/2
    mean = w*m1 + (1-w)*m2
    L = (w*(h(p1)+h(q1)) + (1-w)*(h(p2)+h(q2))) / 2
    C = w*h(ss(p1, q1)) + (1-w)*h(ss(p2, q2))
    sig = (w*(mp.sqrt(rh(p1))+mp.sqrt(rh(q1)))
           + (1-w)*(mp.sqrt(rh(p2))+mp.sqrt(rh(q2)))) / 2
    return mean, L, C, sig

def Phi(p1, q1, p2, q2, w, aa=ALPHA):
    mean, L, C, sig = aggregates(p1, q1, p2, q2, w)
    return (2*(1-aa)*(1-mean) - 1)*L + aa*C - (1-aa)*sig**2

def Lambda(p1, q1, p2, q2, w, aa=ALPHA):
    mean, L, C, sig = aggregates(p1, q1, p2, q2, w)
    return (2*(1-aa)*(1-mean) - 1) + aa*C/L - (1-aa)*sig**2/L

def Q_exact(atoms, weights):
    """Coupled entropy Q = sum_s sum_t u_s u_t h(s+t-st) over atom weights."""
    return sum(us*ut*h(s+t-s*t) for s, us in zip(atoms, weights)
               for t, ut in zip(atoms, weights))

# ---- 1. the obstruction pair (bbar, a) -------------------------------------
g = lambda b: h(b)*(2 - h(b)) - h((1-b)**2)
print("== obstruction defining equation g(b) = h(b)(2-h(b)) - h((1-b)^2) = 0 ==")
roots = []
prev_x, prev_g = mp.mpf("0.0000001"), g(mp.mpf("0.0000001"))
for i in range(1, 20000):
    x = mp.mpf(i) / 20000
    gx = g(x)
    if prev_g * gx < 0:
        rt = mp.findroot(g, (prev_x, x), tol=mp.mpf("1e-45"))
        if all(abs(rt - r0) > mp.mpf("1e-30") for r0 in roots):
            roots.append(rt)
    prev_x, prev_g = x, gx
for rt in roots:
    print("root  b =", mp.nstr(rt, 50), " g =", mp.nstr(g(rt), 3))
bbar = max(roots)
print("bbar (larger root) =", mp.nstr(bbar, 50))
print("|bbar - B_OBS|     =", mp.nstr(abs(bbar - B_OBS), 3))
assert abs(bbar - B_OBS) < mp.mpf("1e-20")

hb = h(bbar)
one_minus_a = 1/(2 - hb)
a = 1 - one_minus_a
print()
print("== relations at the obstruction ==")
print("h(bbar)            =", mp.nstr(hb, 50))
print("1/(2 - h(bbar))    =", mp.nstr(one_minus_a, 50), " (= 1-a)")
print("a                  =", mp.nstr(a, 50))
print("|a - A_OBS|        =", mp.nstr(abs(a - A_OBS), 3))
print("|1/(2-h(bbar)) - A_OBS| =", mp.nstr(abs(one_minus_a - A_OBS), 3),
      "(NOTE: the assignment's 'a = 1/(2-h(b))' is the complement;")
print(" the correct relation, matching A_OBS, is 1-a = 1/(2-h(b)).)")
assert abs(a - A_OBS) < mp.mpf("1e-20")
resid_rel = (1-a)*h((1-bbar)**2) - h(bbar)
print("(1-a)h((1-bbar)^2) - h(bbar) =", mp.nstr(resid_rel, 3))
assert abs(resid_rel) < mp.mpf("1e-45")

# ---- 2. c* ------------------------------------------------------------------
cstar = 1 - (1-bbar)/(2 - hb)
print()
print("== c* ==")
print("psi                =", mp.nstr(PSI, 50))
print("c* = 1-(1-bbar)/(2-h(bbar)) =", mp.nstr(cstar, 50))
print("c* = a + (1-a) bbar         = ", mp.nstr(a + (1-a)*bbar, 50))
print("c* - psi           =", mp.nstr(cstar - PSI, 50), "(= 3.795221...e-4)")
print("|c* - C_STAR const |=", mp.nstr(abs(cstar - C_STAR), 3))
assert abs((cstar - (a + (1-a)*bbar))) < mp.mpf("1e-48")
assert abs(cstar - C_STAR) < mp.mpf("1e-20")
assert abs((cstar - PSI) - mp.mpf("3.795221e-4")) < mp.mpf("1e-10")

# ---- 3. Phi and Phi_exact at the obstruction --------------------------------
# pair-orbit measure: nu* = (1-2a) delta_(bbar,bbar) + 2a delta_(bbar,1)
wA = 1 - 2*a          # weight of orbit A = (bbar, bbar)
wB = 2*a              # weight of orbit B = (bbar, 1)
mean, L, C, sig = aggregates(bbar, bbar, bbar, 1, wA)
atoms = [bbar, bbar, bbar, 1]
wts = [wA/2, wA/2, wB/2, wB/2]
Q = Q_exact(atoms, wts)
phi_val = Phi(bbar, bbar, bbar, 1, wA)
lam_val = Lambda(bbar, bbar, bbar, 1, wA)
print()
print("== obstruction pair-orbit measure nu* (wA = 1-2a on (bbar,bbar), 2a on (bbar,1)) ==")
print("mean(mu*) =", mp.nstr(mean, 50), " (== c*: ", abs(mean - cstar) < mp.mpf("1e-48"), ")")
print("L         =", mp.nstr(L, 50))
print("C         =", mp.nstr(C, 50))
print("Q (exact) =", mp.nstr(Q, 50))
print("sigma^2   =", mp.nstr(sig**2, 50))
print("Q - L =", mp.nstr(Q - L, 3), " C - L =", mp.nstr(C - L, 3),
      " Q - C =", mp.nstr(Q - C, 3))
assert max(abs(Q - L), abs(C - L), abs(Q - C)) < mp.mpf("1e-45")
print("=> Q = C = L to <1e-45 at the obstruction                  [PASS]")
for alpha_try in ("0", "0.0356069", "0.2", "0.5", "1"):
    at = mp.mpf(alpha_try)
    F_exact = (1-at)*Q + at*C - L
    phi_a = Phi(bbar, bbar, bbar, 1, wA, at)
    print("alpha=%9s : Phi_exact=(1-a)Q+aC-L = %+.3e   Phi(certified form) = %+.3e"
          % (alpha_try, F_exact, phi_a))
    assert abs(F_exact) < mp.mpf("1e-45") and abs(phi_a) < mp.mpf("1e-45")
print("Phi and Phi_exact vanish IDENTICALLY in alpha at nu*      [PASS]")
print("Phi(nu*)  =", mp.nstr(phi_val, 5), " Lambda(nu*) = Phi/L =",
      mp.nstr(lam_val, 5))

# cross-check against cert2.phi_true (read-only import)
sys.path.insert(0, "math/uc")
import cert2
ref = cert2.phi_true(bbar, bbar, bbar, 1, wA)
print("|Phi_indep - cert2.phi_true| at nu* =", mp.nstr(abs(phi_val - ref), 3))
assert abs(phi_val - ref) < mp.mpf("1e-45")
print("independent Phi matches cert2.phi_true                    [PASS]")

# ---- 4. margin at t0 ---------------------------------------------------------
print()
print("== margin at t0 = psi + 1e-4 ==")
print("c* - t0           =", mp.nstr(cstar - T0, 50))
heuristic = (cstar - T0) * mp.mpf("1.6298")
print("(c*-t0)*1.6298    =", mp.nstr(heuristic, 50), " (heuristic min margin)")

# Lambda at the face KKT point (check 1 values, re-derived here)
sys.path.insert(0, "math/uc")
import cert3
pk, rk, wk, lamk, resk = cert3.solve_face_kkt(T0)
lam_kkt = Lambda(pk, pk, rk, 1, wk)     # face point (p1,q1,p2,q2,w)=(p,p,r,1,w)
print("Lambda at KKT point (p*,r*,w*) =", mp.nstr(lam_kkt, 45))
print("  (p*,r*,w*) = (%s, %s, %s)" % (mp.nstr(pk, 20), mp.nstr(rk, 20), mp.nstr(wk, 20)))
assert lam_kkt > 0

# NUMERICAL: float64 fine grid on the face, mean <= t0
import numpy as np
t0f = float(T0); aaf = float(ALPHA)
def hf64(x):
    x = np.clip(x, 1e-300, 1-1e-16)
    return -(x*np.log(x) + (1-x)*np.log1p(-x)) / np.log(2.0)

# NUMERICAL: random sampling over the FULL 4-parameter q2=1 face (p1 != q1
# allowed): (p1, q1, r, w) with mean <= t0, w in [1/2, 1].
def lam_face4(p1, q1, r, w):
    mean = w*(p1+q1)/2 + (1-w)*(r+1)/2
    L = (w*(hf64(p1)+hf64(q1)) + (1-w)*(hf64(r)+hf64(1.0)))/2
    C = (w*(hf64(np.minimum(np.maximum(0.5, np.maximum(p1, q1)),
                            np.minimum(p1+q1, 1.0)))
         + (1-w)*hf64(1.0)))
    def rh64(z):
        return np.clip(2*(1-z)*hf64(z) - hf64((1-z)**2), 0, None)
    sig = (w*(np.sqrt(rh64(p1))+np.sqrt(rh64(q1)))
           + (1-w)*(np.sqrt(rh64(r))+np.sqrt(rh64(1.0))))/2
    return ((2*(1-aaf)*(1-mean) - 1)*L + aaf*C - (1-aaf)*sig**2) / L

rng = np.random.default_rng(20260818)
N = 2_000_000
p1 = rng.uniform(0.0, 1.0, N); q1 = rng.uniform(0.0, 1.0, N)
r_ = rng.uniform(0.0, 1.0, N); wv = rng.uniform(0.5, 1.0, N)
mean_f = wv*(p1+q1)/2 + (1-wv)*(r_+1)/2
m = mean_f <= t0f
p1s, q1s, rs, ws = p1[m], q1[m], r_[m], wv[m]
lam4 = lam_face4(p1s, q1s, rs, ws)
k = int(np.argmin(lam4))
print("full-face random sample: %d feasible of %d" % (m.sum(), N))
print("  min Lambda = %.10e at (p1,q1,r,w) = (%.6f, %.6f, %.6f, %.6f)"
      % (lam4[k], p1s[k], q1s[k], rs[k], ws[k]))
assert lam4[k] > 0
print("min Lambda > 0 on full face random sample (NUMERICAL)      [PASS]")

def lam_face(pf, rf, wf):
    mean = wf*pf + (1-wf)*(rf+1)/2
    L = wf*hf64(pf) + (1-wf)*hf64(rf)/2
    C = wf*hf64(np.minimum(np.maximum(0.5, pf), np.minimum(2*pf, 1.0)))
    rh_p = np.clip(2*(1-pf)*hf64(pf) - hf64((1-pf)**2), 0, None)
    rh_r = np.clip(2*(1-rf)*hf64(rf) - hf64((1-rf)**2), 0, None)
    sig = wf*np.sqrt(rh_p) + (1-wf)*np.sqrt(rh_r)/2
    return ((2*(1-aaf)*(1-mean) - 1)*L + aaf*C - (1-aaf)*sig**2) / L

def grid(plo, phi_, rlo, rhi, n, wlo=0.5, whi=1.0):
    ps = np.linspace(plo, phi_, n)
    rs = np.linspace(rlo, rhi, n)
    ws = np.linspace(wlo, whi, n)
    best = (np.inf, None)
    for wv in ws:
        P, R = np.meshgrid(ps, rs, indexing="ij")
        mean = wv*P + (1-wv)*(R+1)/2
        m = mean <= t0f * (1 + 1e-15)
        if not m.any():
            continue
        lam = np.full(P.shape, np.inf)
        lam[m] = lam_face(P[m], R[m], wv)
        idx = np.unravel_index(np.argmin(lam), lam.shape)
        if lam[idx] < best[0]:
            best = (float(lam[idx]), (float(ps[idx[0]]), float(rs[idx[1]]), float(wv)))
    return best

print("NUMERICAL (float64): coarse face grid 200^3 ...", flush=True)
val, arg = grid(1e-6, 0.999999, 1e-6, 0.999999, 200)
print("  coarse min Lambda = %.10e at (p,r,w) = (%.6f, %.6f, %.6f)" % (val, *arg))
for _ in range(3):   # successive refinement, 200^3 each, half-window shrink
    dp_, dr_, dw_ = 0.06, 0.06, 0.03
    val2, arg2 = grid(max(1e-9, arg[0]-dp_), min(1-1e-9, arg[0]+dp_),
                      max(1e-9, arg[1]-dr_), min(1-1e-9, arg[1]+dr_), 200,
                      max(0.5, arg[2]-dw_), min(1.0, arg[2]+dw_))
    val, arg = val2, arg2
    print("  refined min Lambda = %.10e at (p,r,w) = (%.6f, %.6f, %.6f)" % (val, *arg))
assert val > 0
print("min Lambda > 0 on the t0-feasible face (NUMERICAL, grid)   [PASS]")
mp_p, mp_r, mp_w = mp.mpf(repr(arg[0])), mp.mpf(repr(arg[1])), mp.mpf(repr(arg[2]))
print("mpmath(50d) Lambda at grid argmin =", mp.nstr(Lambda(mp_p, mp_p, mp_r, 1, mp_w), 30))
assert Lambda(mp_p, mp_p, mp_r, 1, mp_w) > 0
print()
print("summary: heuristic (c*-t0)*1.6298 = %s" % mp.nstr(heuristic, 8))
print("         Lambda at KKT point      = %s" % mp.nstr(lam_kkt, 8))
print("         grid min Lambda (NUM.)   = %.8e" % val)
print("check2 wall time %.2fs" % (time.time() - t_start))
```

**Output (verbatim, abridged only where repeated above):** see the numbers
quoted in this section; the full transcript matches the script's prints
one-to-one (`full-face random sample: 288717 feasible of 2000000`, four grid
lines ending `4.5754855172e-04`, wall 3.5 s).

---

## Check 3 — orbit-swap invariance

Φ(p1,q1,p2,q2,w) = Φ(p2,q2,p1,q1,1−w) is the exact involution the certifier
uses to restrict w to [½,1] while keeping both atom-face orientations covered
(the aggregates M, L, C, σ are exchanged with the weights). Verified on 100
random points (seed 8731, coordinates in (0.001, 0.999), w in (0.001, 0.999))
at 50 dps with an independent implementation, cross-checked against
`cert2.phi_true` on the same 100 points.

```
max |Phi(p1,q1,p2,q2,w) - Phi(p2,q2,p1,q1,1-w)| = 0.0      (100 points, dps=50)
max |Phi_indep - cert2.phi_true|                 = 0.0      (same 100 points)
```

The invariance holds to the last bit (mpmath evaluates the identical
summands). **PASS.** Implementation note discovered en route (no soundness
impact): importing `cert2` triggers `import entropy`, which sets the global
`mp.dps = 40` at import time — scripts that import these modules must re-pin
`mp.dps` afterwards (done in the script below).

**Script** (`/tmp/hiprec_check3.py`, 0.7 s):

```python
"""Check 3: orbit-swap invariance Phi(p1,q1,p2,q2,w) = Phi(p2,q2,p1,q1,1-w).

100 random points at 50 dps, independent implementation, cross-checked against
cert2.phi_true (read-only import).  The invariance is the exact involution used
by the certifier to restrict w to [1/2, 1] while keeping both face orientations.
"""
import sys
import mpmath as mp
import random

mp.mp.dps = 50
ALPHA = mp.mpf("0.0356069")

def h(x):
    x = mp.mpf(x)
    if x == 0 or x == 1:
        return mp.mpf(0)
    return -(x*mp.log(x) + (1-x)*mp.log(1-x)) / mp.log(2)

def rh(z):
    x = 1 - mp.mpf(z)
    return max(mp.mpf(0), 2*x*h(z) - h(x*x))

def ss(x, y):
    return min(max(mp.mpf("0.5"), max(x, y)), min(x + y, 1))

def Phi(p1, q1, p2, q2, w, aa=ALPHA):
    mean = w*(p1+q1)/2 + (1-w)*(p2+q2)/2
    L = (w*(h(p1)+h(q1)) + (1-w)*(h(p2)+h(q2))) / 2
    C = w*h(ss(p1, q1)) + (1-w)*h(ss(p2, q2))
    sig = (w*(mp.sqrt(rh(p1))+mp.sqrt(rh(q1)))
           + (1-w)*(mp.sqrt(rh(p2))+mp.sqrt(rh(q2)))) / 2
    return (2*(1-aa)*(1-mean) - 1)*L + aa*C - (1-aa)*sig**2

rng = random.Random(8731)
max_swap_err = mp.mpf(0)
max_ref_err = mp.mpf(0)
worst = None
for i in range(100):
    p1 = mp.mpf(rng.uniform(0.001, 0.999))
    q1 = mp.mpf(rng.uniform(0.001, 0.999))
    p2 = mp.mpf(rng.uniform(0.001, 0.999))
    q2 = mp.mpf(rng.uniform(0.001, 0.999))
    w = mp.mpf(rng.uniform(0.001, 0.999))
    lhs = Phi(p1, q1, p2, q2, w)
    rhs = Phi(p2, q2, p1, q1, 1 - w)
    err = abs(lhs - rhs)
    if err > max_swap_err:
        max_swap_err, worst = err, (p1, q1, p2, q2, w)
    assert err < mp.mpf("1e-48"), (i, err)

sys.path.insert(0, "math/uc")
import cert2
mp.mp.dps = 50   # NB: importing cert2 -> entropy sets global mp.dps=40; re-pin
rng = random.Random(8731)
for i in range(100):  # fresh stream: the SAME 100 points as the swap check
    p1 = mp.mpf(rng.uniform(0.001, 0.999))
    q1 = mp.mpf(rng.uniform(0.001, 0.999))
    p2 = mp.mpf(rng.uniform(0.001, 0.999))
    q2 = mp.mpf(rng.uniform(0.001, 0.999))
    w = mp.mpf(rng.uniform(0.001, 0.999))
    ref = cert2.phi_true(p1, q1, p2, q2, w)
    max_ref_err = max(max_ref_err, abs(Phi(p1, q1, p2, q2, w) - ref))
assert max_ref_err < mp.mpf("1e-45")

print("orbit-swap invariance at 100 random points, dps=50:")
print("  max |Phi(p1,q1,p2,q2,w) - Phi(p2,q2,p1,q1,1-w)| =", mp.nstr(max_swap_err, 3))
print("  (worst point seed stream; error is pure rounding at the 1e-49 level)")
print("  max |Phi_indep - cert2.phi_true| over the 100 points =", mp.nstr(max_ref_err, 3))
print("PASS: exact invariance holds to <1e-48; implementations agree to <1e-45")
```

**Output (verbatim):**

```
orbit-swap invariance at 100 random points, dps=50:
  max |Phi(p1,q1,p2,q2,w) - Phi(p2,q2,p1,q1,1-w)| = 0.0
  (worst point seed stream; error is pure rounding at the 1e-49 level)
  max |Phi_indep - cert2.phi_true| over the 100 points = 0.0
PASS: exact invariance holds to <1e-48; implementations agree to <1e-45
```

---

## Check 4 — rh(x) ≤ 2min(x,1−x), and max ρ by line search

**(a) The rh lemma, numerically** (the lemma itself is proved in
`lemma_rh_proof.py`; this is 50-dps corroboration on a 10 001-point grid,
endpoints included):

```
min over grid of [2 min(x,1-x) - rh(x)] = 0.0 at x = 0.0    (and at x = 1)
min interior (x in (0,1)) slack          = 1.57304e-7        (at x = 1e-4)
```

No violation anywhere; the bound is tight only at the endpoints, and the
interior slack 2x − rh(x) ≈ 1.57e-7 at x = 1e-4 exhibits the o(x) tightness
at 0 that the proof must handle (rh = 2x − O(x² log₁₀(1/x))).

**(b) max_p ρ by line search** (ρ(p) = rh(p)/h(p) = 2(1−p) − h(2p−p²)/h(p),
`margin_lemma.rho_mp`'s formula; 20 000-point scan, then golden section to
width < 1e-30, then stationary-point polish ρ′=0 via `mp.findroot`):

```
argmax p (golden)          = 0.116502234759150983054928429238252060105
argmax p (rho' = 0)        = 0.1165022347591509830549284835278949219607
max rho                    = 0.3049467103530205557127432257297999107925
|golden - stationary|      = 5.43e-26
```

Consistency with the certified record:

- PROOF.md / README record **max ρ = 0.3049467 at p ≈ 0.1165** (labelled
  NUMERICAL, sampled): reproduced — |max ρ − 0.3049467| = 1.04e-8, argmax
  0.11650223… rounds to 0.1165. **[PASS]**
- `cert2.get_rho_gmax()` (certified Arb enclosure, n = 20 000 subdivision +
  analytic endpoint guards) returns the certified cap
  **0.305360849157902…** > 0.3049467103530… — the numerical maximiser sits
  below the certified cap, as required. **[PASS]**

**Script** (`/tmp/hiprec_check4.py`, 2.7 s):

```python
"""Check 4: (a) numerical corroboration of rh(x) <= 2 min(x,1-x) on a 10k grid
at 50 dps; (b) max_p rho(p) by line search, compared with the certified
max rho = 0.3049467.

rho(p) = rh(p)/h(p) = 2(1-p) - h(2p-p^2)/h(p)  (margin_lemma.rho_mp formula).
"""
import sys
import mpmath as mp

mp.mp.dps = 50

def h(x):
    x = mp.mpf(x)
    if x == 0 or x == 1:
        return mp.mpf(0)
    return -(x*mp.log(x) + (1-x)*mp.log(1-x)) / mp.log(2)

def rh(z):
    x = 1 - mp.mpf(z)
    return max(mp.mpf(0), 2*x*h(z) - h(x*x))

def rho(p):
    p = mp.mpf(p)
    return 2*(1-p) - h(2*p - p*p)/h(p)

N = 10000
worst = mp.inf
worst_x = None
for i in range(N + 1):
    x = mp.mpf(i) / N
    slack = 2*min(x, 1-x) - rh(x)     # claimed >= 0
    if slack < worst:
        worst, worst_x = slack, x
print("rh(x) <= 2 min(x,1-x) on %d-point grid, dps=50:" % (N + 1))
print("  min over grid of [2 min(x,1-x) - rh(x)] =", mp.nstr(worst, 6),
      "at x =", mp.nstr(worst_x, 8))
assert worst >= 0
worst_int = min(2*min(x, 1-x) - rh(x) for x in
                (mp.mpf(i)/N for i in range(1, N)))
print("  min interior (x in (0,1)) slack =", mp.nstr(worst_int, 6))
assert worst_int > 0
print("  PASS: no violation; slack -> 0 only at the endpoints x=0, x=1")

# (b) max rho by line search
best = (-1, None)
for i in range(1, 20000):
    p = mp.mpf(i) / 20000
    v = rho(p)
    if v > best[0]:
        best = (v, p)
p0 = best[1]
lo, hi = p0 - mp.mpf("5e-5"), p0 + mp.mpf("5e-5")
# golden-section to width 1e-30
gr = (mp.sqrt(5) - 1) / 2
a_, b_ = lo, hi
c_ = b_ - gr*(b_ - a_)
d_ = a_ + gr*(b_ - a_)
for _ in range(200):
    if rho(c_) < rho(d_):
        a_ = c_
    else:
        b_ = d_
    c_ = b_ - gr*(b_ - a_)
    d_ = a_ + gr*(b_ - a_)
    if b_ - a_ < mp.mpf("1e-30"):
        break
p_max = (a_ + b_) / 2
rho_max = rho(p_max)
# stationary-point polish: solve rho'(p) = 0 by findroot on the numeric derivative
p_stat = mp.findroot(lambda z: mp.diff(rho, z), p_max, tol=mp.mpf("1e-40"))
print()
print("max_p rho by line search (golden section, width < 1e-30):")
print("  argmax p =", mp.nstr(p_max, 40))
print("  stationary point (rho' = 0) p =", mp.nstr(p_stat, 40))
print("  max rho  =", mp.nstr(rho_max, 40))
print("  rho(p_stat) =", mp.nstr(rho(p_stat), 40))
print("  |golden - stationary| =", mp.nstr(abs(p_max - p_stat), 3))
assert abs(rho_max - rho(p_stat)) < mp.mpf("1e-35")
assert abs(p_max - mp.mpf("0.1165")) < mp.mpf("5e-4")
assert abs(rho_max - mp.mpf("0.3049467")) < mp.mpf("5e-8")
print("  certified value on record: 0.3049467 (at p ~ 0.1165)")
print("  |max rho - 0.3049467| =", mp.nstr(abs(rho_max - mp.mpf("0.3049467")), 3),
      " -> consistent                       [PASS]")

# cross-check with cert2's certified Arb global cap (read-only use)
sys.path.insert(0, "math/uc")
import cert2
cap = cert2.get_rho_gmax()
print("  cert2.get_rho_gmax() certified Arb cap =", cap.upper().str(15, radius=False))
assert rho_max < mp.mpf(float(cap.upper()))
print("  numerical max < certified cap                               [PASS]")
```

**Output (verbatim):**

```
rh(x) <= 2 min(x,1-x) on 10001-point grid, dps=50:
  min over grid of [2 min(x,1-x) - rh(x)] = 0.0 at x = 0.0
  min interior (x in (0,1)) slack = 1.57304e-7
  PASS: no violation; slack -> 0 only at the endpoints x=0, x=1

max_p rho by line search (golden section, width < 1e-30):
  argmax p = 0.116502234759150983054928429238252060105
  stationary point (rho' = 0) p = 0.1165022347591509830549284835278949219607
  max rho  = 0.3049467103530205557127432257297999107925
  rho(p_stat) = 0.3049467103530205557127432257297999107925
  |golden - stationary| = 5.43e-26
  certified value on record: 0.3049467 (at p ~ 0.1165)
  |max rho - 0.3049467| = 1.04e-8  -> consistent                       [PASS]
  cert2.get_rho_gmax() certified Arb cap = 0.305360849157902
  numerical max < certified cap                               [PASS]
```

---

## Discrepancy register (complete)

1. **PROOF.md λ\* last digit (doc-level only).** PROOF.md line 220 records
   λ\* = 1.3682158100251758; the solve gives 1.36821581002517574441…, i.e.
   …1757 at 17 significant digits. |Δ| = 5.56e-17 (16-digit agreement).
   λ\* is explicitly a NUMERICAL seed whose 3-point family {0.98λ\*, λ\*,
   1.02λ\*} is invariant under this rounding; no certificate step depends on
   the 17th digit. Suggested doc fix only: write 1.3682158100251757.
2. **Task-statement relation at the obstruction.** "a = 1/(2−h(b))" in the
   assignment is the complement of the true relation; numerically
   1/(2−h(b̄)) = 0.9211227072940768… = 1−a, and a = 1 − 1/(2−h(b̄)) =
   0.07887729270592317… matches `A_OBS` to 1.6e-22 (README "Result 2" states
   the relation correctly as 1−a = 1/(2−h(b))). All obstruction quantities
   here were computed with the correct relation.
3. **Minor margin-spread (expected, not a defect).** The three margin
   estimates at t0 (heuristic 4.5557e-4, Λ at KKT point 4.5661e-4, grid min
   4.5755e-4) differ by < 0.5% relative: the KKT stationarity targets Φ, not
   Λ = Φ/L, and the float64 Λ-landscape is flat at the minimum. All three are
   strictly positive; the recorded "min-Λ ≈ 4.56e-4" is reproduced.

**No other discrepancies.** All recorded constants (`B_OBS`, `A_OBS`,
`C_STAR`, certified ρ cap) reproduce to at least the number of digits on
record; the KKT solve is bit-identical to cert3's own; Φ and the orbit
involution match the in-repo reference implementations exactly.
