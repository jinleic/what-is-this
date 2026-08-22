# The paired-momentum defect: the exact weakest additional inequality class for the upper endpoint

**File classification.** [THEOREM] Theorem PM (the gap theorem): for every
`K* < I3/2` the *paired-momentum ceiling-defect floor* —
`liminf_{L even} Delta_L(K*) > delta*(K*) := I3/(2K*) - 1` — is exactly the
weakest additional inequality that certifies `K_c <= K*`; it is equivalent,
with **no** finite-size correction, to a uniform magnetisation floor, and no
consequence of the audited two-point constraint system can supply it
(`proofs/upper_infrared.md`, Theorem D).  [COMPUTATION] exact instance data
on every torus the in-house engine reaches (`2x2x2`, `2x2x3`, `2x2x4`,
`3x3x2`), certified massive-Watson thresholds `I3(mu^2)` and
Brillouin-zone shell thresholds `eps*(K*, Lambda)`.

Producer: `experiments/e145_paired_momentum.py`.  Standalone auditor:
`tests/test_paired_momentum.py` (full-configuration enumeration, reciprocal
exponential route, binomial-product walk counts, own interval convolution,
own shell grid).  Artifact: `results/bounds/paired_momentum.json`.

**Benchmark discipline.** `K_c = 0.221654626` was not used to select any
threshold, target, or rounding.  Evaluation points are the certified
endpoints (the `I3/2` bracket, the certified lower floor from
`results/bounds/saw_union4.json`) and the plain rationals `63/250`,
`1/4`, `49/200`, `6/25`, `11/50`.

---

## 0. Context: the two walls this front starts from

1. **Peierls wall** (`proofs/upper_beyond.md`, Theorem P): the plain
   contour-sum certificate fails for every `K <= K*` with certified
   `K* = 0.25583... > I3/2` — margin `V(x_-) - 1/2 = 0.08832` at the
   incumbent, tail mass dominating any head gain.  The head-extension
   programme is closed.
2. **Two-point wall** (`proofs/upper_infrared.md`, Theorem D): the
   reflection-positivity/GKS/infrared constraint system `S(K)` is exactly
   saturated at `2K = I3` by `h_*(k) = 1/(I3*lambda(k))` with zero
   long-range order, and the finite-torus LP floors evaporate as
   `5168/(525L)` (`proofs/mag_floor.md`).  The two-point class is closed.

The assignment: identify the **weakest additional inequality class** that
could break `I3/2`, state it as an exact sufficient condition, and test the
simplest instance on exact small-lattice data.

## 1. The master identity and the gap theorem

Let `T_L` be any finite torus with site set of size `N`, `h = 0`,
`M = sum sigma`, `M_L^2 = <M^2>/N^2`, `G_L(z) = <sigma_0 sigma_z>`, and
`Ghat_L(k) = sum_z G_L(z) cos(k.z)` its (real) Fourier transform; for a
torus with bond multiplicity `mult_i` in direction `i` (a periodic
direction of length 2 carries two parallel bonds, `notes/CODE_API.md`)
the graph dispersion is `lambda_T(k) = sum_i mult_i (1 - cos k_i)`.

**[LEMMA M1, exact identity]** On any finite torus,

```
  1 = M_L^2 + (1/N) sum_{k != 0} Ghat_L(k)                       (Parseval)
  Delta_L(K) := (1/N) sum_{k != 0} [ 1/(2K lambda_T(k)) - Ghat_L(k) ]
             = M_L^2(K) - 1 + C_L(0)/(2K),
  C_L(0) := (1/N) sum_{k != 0} 1/lambda_T(k).
```

*Proof.* `sigma_x^2 = 1` and the mode orthogonality give Parseval; the
second line is Parseval plus termwise algebra — no inequality is used, and
no reflection-positivity hypothesis enters.  On even tori the
Fröhlich–Simon–Spencer ceiling makes every bracket nonnegative, hence
`Delta_L >= 0`; on the small tori here this is a *measured* statement, not
a theorem (§3).  `square`

`Delta_L` is the **paired-momentum ceiling defect**: the amount by which
the two-point spectral profile sits below the infrared ceiling, mode by
mode.  The `I3`-anchored version

```
  DeltaHat_L(K*) := M_L^2(K*) - 1 + I3/(2K*)
```

differs from `Delta_L` by `(C_L(0) - I3)/(2K*) = o(1)` (the even-torus
Riemann sums converge to `I3`; lower bound `C_L(0) >= I3 - A/L` with
`A = 992015/408608` is certified in `proofs/mag_floor.md` §2.3).

> **[THEOREM PM — the gap theorem]**  Fix `0 < K* < I3/2` and `eta > 0`,
> and let `delta*(K*) := I3/(2K*) - 1 = (I3 - 2K*)/(2K*) > 0`.
> The following are equivalent for even tori:
>
> 1. `DeltaHat_L(K*) >= delta*(K*) + eta` for all even `L >= L0`;
> 2. `M_L^2(K*) >= eta` for all even `L >= L0`.
>
> Either implies `K_c <= K*`.  Moreover `delta*(K*)` is the **exact**
> threshold: with `liminf DeltaHat_L(K*) = delta*(K*)` precisely, only
> `M_L^2 >= 0` follows, and no certificate is obtained.

*Proof.* `M_L^2 = DeltaHat_L + 1 - I3/(2K*)` is an identity, and
`1 - I3/(2K*) + delta*(K*) = 0` exactly; so (1) and (2) are the same
statement.  A uniform positive lower bound on `M_L^2` at `K*` forces
positive spontaneous magnetisation at `K*` by the standard finite-volume
argument *inherited verbatim from the audited infrared route*
(`proofs/upper_infrared.md`, eq. (4)): `M_L^2 = (1/N) sum_z G_L(z) <=
max_z G_L(z)`, a pigeonhole over the `(2R+1)^3` offsets of radius `R`
forces arbitrarily distant offsets with `G_per(z) >= eta/2`, boundary
monotonicity (`plus >= periodic >= free`, Griffiths) transfers them to the
plus state, and hence `K* >= K_c`.  The
threshold is sharp because the barrier profile of Theorem D has
`DeltaHat equiv delta*`-convergent values with `M^2 = 0`: it is feasible
for the two-point system at exactly this defect level.  `square`

**Weakest-class statement.** Any inequality — of whatever proof origin —
that improves the upper endpoint *is* a lower bound on
`liminf DeltaHat_L(K*)` above `delta*(K*)`: by the identity, an endpoint
improvement at `K*` means `M_L^2(K*) >= eta > 0`, which *is* the defect
floor.  The search for the weakest additional inequality class therefore
has a single exact target:

> **[CONDITION Delta-FLOOR]**  *If there exist `eta > 0` and `L0` with*
> `DeltaHat_L(K*) = M_L^2(K*) - 1 + I3/(2K*) >= (I3-2K*)/(2K*) + eta`
> *for all even `L >= L0`, then `K_c <= K*`.*

### 1.1 The exact threshold table `[COMPUTATION, exact rational intervals]`

With the certified `I3` interval of `results/bounds/upper_infrared.json`:

| `K*` | `delta*(K*) = I3/(2K*)-1` |
|---|---|
| `63/250 = 0.252` (just below the incumbent `0.2527310098...`) | `0.0029008...` |
| `1/4 = 0.25` | `0.0109240...` |
| `49/200 = 0.245` | `0.0315551...` |
| `6/25 = 0.24` | `0.0530459...` |
| `11/50 = 0.22` | `0.1487773...` |
| certified lower floor `0.2122159753...` | `0.1909142...` |

So: to beat the incumbent by `0.00073` one needs a uniform defect floor of
`0.0029`; by `0.00273` (`K* = 1/4`), `0.01092`; to reach the certified
*lower* endpoint, `0.19091`.  At the incumbent itself the threshold is
`0` — which is exactly why the two-point class, whose barrier realises
defect `o(1)`, cannot move.

## 2. The available inequality classes, each with its exact condition

### (A) Two-point GKS/LP class — exhausted

Exact sufficient condition of this class: `exists eta>0, L0` with
`F(K*, L) = 1 - S_*(K*, L) >= eta` for all even `L >= L0` (`S_*` the
finite-torus Parseval/infrared/GKS/cap/energy LP of `proofs/mag_floor.md`).
**Refuted** for every `K <= I3/2`: `F(K,L) <= 5168/(525L)`
(`proofs/mag_floor.md`, Theorem (3)).  Theorem D excludes every
rearrangement.  Nothing remains here; all further classes are strictly
four-point or field-derivative objects.

### (B1) Massive-IR (mode-wise four-point-derivable condition)

> **[CONDITION PM-mu]** *If there exist `mu > 0` and `L0` such that for all
> even `L >= L0` and all `k != 0`:*
> `Ghat_L(k) <= 1/(2K* (lambda(k) + mu^2))`, *and* `I3(mu^2) < 2K*`,
> *then `K_c <= I3(mu^2)/2 <= K*`.*

*Proof.* Sum the ceiling in Parseval; the `k=0` exclusion costs
`(1/N)/(2K* mu^2) -> 0`, and
`(1/N) sum_{k!=0} 1/(lambda+mu^2) -> I3(mu^2)`.  `square`

Certified `[COMPUTATION]` (directed positive interval series over the exact
walk return probabilities `R_{2n} = N_{2n}/6^{2n}`,
`I3(mu^2) = sum_n 3^n R_n / (3+mu^2)^{n+1}`, tail `<= r^{N+1}/mu^2` since
`R_n <= 1`):

| `mu^2` | certified `I3(mu^2)` | `I3(mu^2)/2` |
|---|---|---|
| `1/2` | `0.421260237549...` | `0.210630118774...` |
| `1/4` | `0.440435872998...` | `0.220217936499...` |
| `1/8` | `0.456571449668...` | `0.228285724834...` |
| `1/16` | `0.469389951217...` | `0.234694975608...` |
| `1/32` | `0.479192407559...` | `0.239596203779...` |
| `1/128` | `>= 0.491862184635` (lower bound only at `N=3200`) | — |
| `1/1024` | `>= 0.500317453224` (lower bound only) | — |

Two certified consequences:

* **Mu-bracket.** `I3(1/1024) > 1/2 > I3(1/32)` gives
  `mu*(1/4)^2 in (1/1024, 1/32)` — i.e. `mu* in (1/32, 1/(4*sqrt2))`:
  a *tiny* uniform mass already certifies `K_c <= 1/4`.  For
  `K* in {49/200, 6/25}`: `mu*^2 in (1/128, 1/32)`.  (Finer brackets need
  the `N ~ 4*10^4` convolution; open this wave.)
* **Self-limitation.** `I3(1/2)/2 = 0.2106301...` is *below* the certified
  lower endpoint `K_c >= 0.2122159753...`, so PM-mu with `mu^2 >= 1/2`
  cannot hold anywhere below `K_c`: any provable massive ceiling near
  criticality must have `mu^2 < 1/2`.  The whole game lives in
  `mu^2 in (0, 1/2)`, tightening toward `(1/1024, 1/32)` as the target
  approaches the incumbent.

### (B2) Shell-relative defect (the mildest quantitative form)

> **[CONDITION PM-S]** *If there exist `eps > 0`, `Lambda > 0`, `L0` such
> that for all even `L >= L0` and all modes with `lambda(k) >= Lambda`:*
> `Ghat_L(k) <= (1 - eps)/(2K* lambda(k))`, *then*
> `DeltaHat_L(K*) >= eps*J(Lambda)/(2K*) + o(1)`, *and the condition*
> `eps > eps*(K*, Lambda) := (I3 - 2K*)/J_lo(Lambda)` *certifies
> `K_c <= K*`.*

*Proof.* The defect on the shell sums to
`(eps/(2K*)) * <1/lambda ; lambda >= Lambda>`, and the shell weight is
lower-bounded by the certified grid value `J_lo(Lambda) <= J(Lambda)`;
insert into Theorem PM.  `square`

`J_lo(Lambda)` is certified by exact rational Brillouin-zone cell sums
(grid `96^3`, monotone cos bounds; the same grid sum gives the machine
check `0.481480... <= I3` for the enclosure machinery):
`J_lo(2) = 0.2385023...`, `J_lo(3) = 0.1241778...`, `J_lo(4) = 0.0438377...`.

Exact thresholds `[COMPUTATION]`:

| `K*` | `eps*(K*, Lambda=2)` | `eps*(K*, Lambda=3)` |
|---|---|---|
| `63/250` | `0.00613...` | `0.01177...` |
| `1/4` | `0.02290...` | `0.04399...` |
| `49/200` | `0.06483...` | `0.12452...` |
| `6/25` | `0.10676...` | `0.20504...` |

A **2.3% uniform relative ceiling defect on all modes with
`lambda(k) >= 2`** certifies `K_c <= 1/4`; a **0.61%** defect certifies
`K_c <= 0.252`.  These are the mildest quantitative targets on record for
this endpoint.

### (B3) What four-point *sign* structure cannot do `[PROPOSITION W]`

No inequality satisfied by every centred Gaussian field with
nonnegative-definite covariance — GKS-I nonnegativity `G(z) >= 0`, the
Lebowitz four-point upper bound (Wick sum), Jensen
`<X^4> >= <X^2>^2`, positive-definiteness — can imply the Delta-FLOOR: the
Gaussian field with spectral density `1/(I3 lambda)` satisfies all of them
(the Lebowitz bound with equality), has `M^2 = 0` and hence defect exactly
at the threshold.  `square` (Immediate from Theorem D's witness; recorded
because it rules out an entire natural class of proposed "four-point
improvements".)  Any useful four-point input must be *violated* by the
Gaussian — i.e. use `sigma_x^2 = 1` (bounded spins) or field-derivative
structure.  The paired-momentum intensities
`<|sigmahat(k)|^4>` stored in the artifact are exactly this kind of data
(§3.3).

### (C) Random-current class (the "Bryn" slot)

No source named "Bryn" exists in `sources/` or in any repo note (searched;
no match), so no audited text is available.  The class is therefore stated
as the standard random-current/switching template with its exact
in-house-integrable condition:

> **[CONDITION PM-ABF(kappa)]** *If for the infinite-volume susceptibility*
> `chi(K) = sum_z G_K(z)` *one has* `d chi/dK >= kappa chi^2` *a.e. on*
> `[K0, K*)` *for some* `kappa > 0`, *and* `chi(K0) >= 1/(kappa (K* - K0))`,
> *then `chi(K*) = infinity`, hence `K_c <= K*`.*

*Proof.* `d(1/chi)/dK <= -kappa` gives
`1/chi(K*) <= 1/chi(K0) - kappa(K*-K0) <= 0`.  `square`
[TheOREM: integration only; the differential inequality itself is
**[EXTERNAL/OPEN-in-house]** — deriving it needs the switching lemma,
not attempted this wave; `kappa` has no certified positive value in-house.
This is the precise open slot of the class.]

### (D) Griffiths–Hurst–Sherman class

GHS (concavity of `m_L(h)` in `h >= 0`) is an [EXTERNAL] theorem for
Ising ferromagnets (no audited full text in the manifest; not re-derived
here).  Its h-direction floors are provably too weak:

> **[PROPOSITION GHS-h, conditional on GHS]** *If `m_L(h0) >= m0 > 0` for
> all `L` at some fixed `h0 > 0`, then only*
> `M_L^2 = chi_L(0)/N >= m0/(h0 N) -> 0` *follows: no uniform floor, no
> endpoint improvement.*

*Proof.* Concavity with `m_L(0) = 0` gives `m_L(h0) <= h0 chi_L(0)`.  `square`

The surviving GHS-class content is derivative (third-cumulant /
`dchi/dh <= 0`) and couples to the same `dchi/dK` structure as (C); GHS
alone supplies no conversion of four-point sign data into a Delta-FLOOR
(PROPOSITION W applies to its Gaussian-saturated instances).

## 3. The exact finite-lattice instance

### 3.1 Engine and anchors `[COMPUTATION]`

Every torus the in-house transfer engine reaches is enumerated exactly
(`sigma_0 = +1` half space, integer signed histograms in `x = e^{-2K}`,
length-2 directions doubled per `notes/CODE_API.md`); `2^17` half-states
for the largest, `3x3x2` (`N = 18`, 54 bonds).  Anchors, all exact:
broken-bond polynomials equal to `ising.transfer_matrix.torus_broken_bond_poly`
for all four tori; Parseval `(1/N) sum_k Ghat = 1`; the paired identity
`N sum_z G(z)^2 = sum_k Ghat(k)^2` on every torus; offset-class orbit
constancy at the signed-integer-histogram level; and the wave-8 audit
anchor `G(e) = 4531/15844`, `Ghat(pi,0,0) = 280125/269348` on the
single-bond 8-site torus at `x = 3/5` reproduced bit-exactly
(`proofs/upper_infrared.md` §3.2).

Evaluation points: the certified incumbent bracket (rational `x*`
endpoints of `results/bounds/upper_beyond.json`, truncated to 60 decimals
*in the containment direction*) and own exact enclosures of `e^{-2K*}` at
`K* in {1/4, 49/200, 6/25, lower floor}`.  All bracket arithmetic is exact
`Fraction` work; transcendental enclosures are alternating-Taylor
(producer) vs. reciprocal-positive-series (auditor).

### 3.2 The defect and the ceiling at the incumbent bracket

`Delta_L` (exact identity of Lemma M1) at the incumbent `K = I3/2`
(intervals of width `10^-60`, shown to 6 places):

| torus | `N` | `Delta_L` | `M_L^2` | modes violating the ceiling |
|---|---|---|---|---|
| `2x2x2` | 8 | `-0.023029` | `0.678152` | 7 of 7 |
| `2x2x3` | 12 | `+0.091985` | `0.614512` | 5 of 11 |
| `2x2x4` | 16 | `+0.106961` | `0.540578` | 10 of 15 |
| `3x3x2` | 18 | `+0.203061` | `0.593201` | **0 of 17** |

Two exact structural facts (per-mode data in the artifact):

1. **Every violating mode on every torus has a nonzero momentum component
   along a length-2 (doubled-bond) direction.**  The degenerate directions
   are exactly where reflection positivity degenerates (the wave-8 audit
   found the same on the single-bond `L=2` torus at `v = 1/4`); the modes
   of non-degenerate directions obey the ceiling even at `N = 12`.
2. On the largest reachable torus `3x3x2` the mode-wise ceiling holds
   strictly at the incumbent on **all 17 nonzero modes**; the smallest
   relative defect is on the degenerate-direction mode `(0,0,1)`
   (`lambda_T = 4`, `Ghat = 0.492027`, defect `+0.00257`), while the
   low-dispersion modes carry huge defects (`lambda_T = 3/2`:
   `Ghat = 0.735511`, defect `+0.583414`, i.e. a 44% relative defect).

Along the size tower `2x2x2 -> 2x2x3 -> 2x2x4 -> 3x3x2` the incumbent
defect rises `-0.023 -> +0.092 -> +0.107 -> +0.203`, an order of magnitude
above `delta*(1/4) = 0.0109` already at `N = 12`.

### 3.3 Four-point / second-moment instance data

The summed two-point second moment `S2 = sum_z G(z)^2` and Binder ratios
`B = <M^4>/<M^2>^2` at the incumbent:

| torus | `S2` | `B` (Gaussian 3, Lebowitz `<= 3`) |
|---|---|---|
| `2x2x2` | `3.812639` | `1.3014` |
| `2x2x3` | `4.726885` | `1.3408` |
| `2x2x4` | `4.987764` | `1.4272` |
| `3x3x2` | `6.540624` | `1.3108` |

and the paired-momentum four-point intensities on the integer-coefficient
modes (e.g. `3x3x2`, mode `(0,0,1)`: `<|sigmahat|^4> = 8.8565...` against
the Gaussian reference `235.31...`, i.e. 3.8% of Gaussian — the mode
intensities of the ±1 field fluctuate far below Gaussian).

> **[PROPOSITION S2 — second-moment floors are vacuous]**  No condition of
> the form `sum_{z in T_L} G_L(z)^2 >= c` (fixed `c`, all large `L`) can
> break `I3/2`: the barrier profiles with spectral density
> `1/(I3 lambda(k))` on `k != 0` have `M^2 = 0` and
> `sum_z G^2 = (1/N) sum_{k!=0} 1/(I3 lambda)^2 >= L/(576 pi^4 I3^2) -> infinity`
> by the dyadic shell count (modes with `max_i |k_i| in (2^{-j-1}pi, 2^{-j}pi]`
> number at least `(L 2^{-j-2})^3` for `L 2^{-j} >= 8`, each contributing
> `(1/N) * 8/(9 pi^4) * 2^{4j}`; summing `j <= log2(L/8)` gives the
> constant, with `pi <= 355/113` exact).  `square`

This closes the literal suggestion "if `sum_x m(x)^2 >= c`" from the
assignment: the correct second-moment-type quantity is the *defect*
`DeltaHat` of Theorem PM, not the raw second moment.

### 3.4 Verdict on the simplest instance

**INDECISIVE — supports plausibility, decides nothing.**  Precisely:

* *Supports:* on the largest reachable torus the mode-wise ceiling holds
  on every non-degenerate mode and on all modes at the incumbent, the
  defect `Delta = 0.2031` exceeds `delta*(1/4)` by a factor 18 and
  `delta*(63/250)` by a factor 70, and the defect grows monotonically
  along the reachable tower.
* *Contradicts (the naive finite-`L` forms):* the mode-wise ceiling and
  `Delta_L >= 0` are **not** finite-`L` facts — they fail on the
  degenerate-direction modes of three of the four tori.  PM-mu and PM-S
  are asymptotic statements; no reachable lattice can confirm or refute a
  `liminf` over even `L`.
* The exact thresholds (§1.1, §2) are unconditional; the open work is
  *deriving* a uniform defect floor from four-point/random-current/GHS
  input, not measuring it.

### 3.5 Why no in-house series order suffices at the incumbent

The HT/LT series artifacts (`results/series/`, free energy to `v^28`, LT to
`x^34`) converge only up to their radius, `v_c = tanh K_c approx 0.218 <
v* = tanh(I3/2) approx 0.2475`.  A series evaluation of `chi`, `S2`, or any
susceptibility-type quantity *at the incumbent bracket* would require
proving divergence at `v*` — which is exactly the open problem
(`chi(v*) = infinity iff K_c <= I3/2`, the endpoint itself).  Partial sums
with nonnegative coefficients are valid lower bounds only inside the disc
of validity, so no available order gives certified instance values at the
incumbent: **exact finite-lattice enumeration is the only in-house route**,
and that is what §3 uses (with the bracket-endpoint error accounting
described there).  No benchmark value enters this reasoning; it uses only
the certified inequality `K_c >= 0.2122159753...` (which already forces
`v_c >= 0.2101`, and the actual obstruction is the radius, not its exact
location).

## 4. What is not claimed; walls reached

* No improved upper endpoint.  Theorem PM is a gap theorem: it identifies
  the exact target and thresholds, and the instance data is finite-volume
  only.
* `C_L(0) -> I3` is inherited standard input (used by `e72` eq. (4)); the
  lower bound is certified in `proofs/mag_floor.md`; the upper side is the
  elementary cell comparison (§1), not re-proved line-by-line here.
* The ABF differential inequality of class (C) has no in-house derivation
  and no audited source text in the manifest; the "Bryn" attribution could
  not be identified anywhere in the repository.  Both are recorded as open
  slots, not as available input.
* The massive-Watson values at `mu^2 <= 1/128` are certified lower bounds
  only at `N = 3200` (the tail bound `r^{N+1}/mu^2` is loose there);
  tightening needs the `N ~ 4*10^4` convolution.
* Engineering note (repo-wide relevance): in this venv (CPython 3.14,
  mpmath 1.3.0) the idiom `mp.dps = n` is **inert** — precision silently
  stays at 53 bits; `mp.mp.dps = n` / `mp.iv.dps = n` must be used.  The
  producer sets both explicitly; scripts relying on `mp.dps` elsewhere
  deserve re-auditing.

## 5. Reproduction

```
PYTHONPATH=src .venv/bin/python experiments/e145_paired_momentum.py
PYTHONPATH=src .venv/bin/python tests/test_paired_momentum.py
```

The auditor recomputes, by independent engines: the full `2^N`
enumerations and all stored finite-lattice values at the incumbent
bracket endpoint (exact `Fraction` equality with a bracket endpoint), the
`delta*` table, the walk counts by the binomial-product identity, its own
interval convolution for `I3(mu^2)` at `mu^2 in {1/2, 1/8}` (overlap
asserted), the grid-48 shell weights with its own cos enclosures, the
massive-threshold bracket logic, the frozen transfer-matrix polynomials,
and the wave-8 anchor.
