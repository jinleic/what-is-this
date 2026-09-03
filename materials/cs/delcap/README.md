# `delcap/` — rigorous enclosure of the binary deletion channel capacity bounds

**Status: BENCHMARK (no campaign run yet).** The capacity $\mathcal C(\mathrm{BDC}_d)$
of the i.i.d. binary deletion channel is open. Every published numeric bound
constant on it comes from a **floating-point, computer-aided Blahut–Arimoto
analysis**, and **none has ever been recomputed in rigorous interval
arithmetic**. This target builds the first certified enclosure table.

## The published chain

Owner-read (**V-owner**, 2026-08-29) — arXiv:2305.07156, Ittai Rubinstein &
Roni Con, *"Improved Upper and Lower Bounds on the Capacity of the Binary
Deletion Channel"* (2023-05-11) states the landscape verbatim:

* lower $\mathcal C(\mathrm{BDC}_d)>0.1185\,(1-d)$ — Mitzenmacher–Drinea
  (ITIT 2006), Kirsch–Drinea (ITIT 2009), via run-length distributions;
* upper $\mathcal C(\mathrm{BDC}_d)<0.4143\,(1-d)$ for $d>0.65$ — Fertonani–Duman
  (ITIT 2010), Dalai (ISIT 2011), Rahmati–Duman (ITIT 2014), *"computer aided
  analyses based on the Blahut-Arimoto algorithm"*;
* **this paper**: upper $\mathcal C(\mathrm{BDC}_d)<0.3745\,(1-d)$ for all
  $d\ge0.68$, and lower $\mathcal C(\mathrm{BDC}_d)>0.1221\,(1-d)$, both by
  extending Blahut–Arimoto to lower space complexity.

Scout-verified (**V-scout**, owner read pending — see pre-campaign):

* arXiv:**0810.0785**, Fertonani & Duman, *"Novel Bounds on the Capacity of the
  Binary Deletion Channel"* (2008-10-04) — four genie-aided upper bounds with
  an infinite series reduced to closed form.
* arXiv:**0912.5176** (2009-12-29) and arXiv:**1104.5546** (2011-04-29),
  Kanoria & Montanari — small-$d$ capacity expansion, two and three leading
  terms; coefficients published as rationals but hand-verified.
* arXiv:**2604.05867**, Pinto–Ribeiro (2026) — $\mathcal C\le0.3578\,(1-d)$ for
  $d\ge0.64$ via a **CUDA/parallelized** Blahut–Arimoto; no code release
  located, so the constant rests on an unreproduced GPU pipeline.
* arXiv:**2504.20961**, Morozov & Duman (2025-04-29) — *"Simple Finite-Length
  Achievability and Converse Bounds for the Deletion Channel and the Insertion
  Channel"*; no independent reproduction found.
* arXiv:**2607.19559**, Tavakoli, Nguyen & Bose (2026-07-21) — q-ary
  combinatorial capacity bounds, with the sandwich derived from *"numerical
  experiments at $n=3,5,10$ and $q=2,3$"* only.

## Why this is a certificate target

Blahut–Arimoto is an alternating maximization over finite probability
simplices; the published constants are the numerical optima of finite programs.
Replacing the inner maximization by a **rigorous interval upper bound**
(Arb balls plus explicit branch-and-bound with a stated Lipschitz/monotonicity
argument) converts each published float into a certified enclosure. The
deliverable is a table of (published constant, certified interval,
PASS/FAIL/IMPROVE) with a reproducible artifact per row. Every outcome is
publishable: uniform confirmation is the first rigorous enclosure table for a
named open capacity; any non-enclosure is a correction to the literature.

## Gates

1. **Gate A (enclosure of the two crux constants — refute-or-confirm, hours).**
   Re-derive Rubinstein–Con's upper $0.3745\,(1-d)$ at $d\in\{0.68,0.75,0.85,0.95\}$
   and lower $0.1221\,(1-d)$, in `python-flint`/Arb with outward rounding, the
   inner Blahut–Arimoto maximization replaced by a certified interval bound.
   **Pass:** each published constant is enclosed within its published rounding
   precision. **Refutation:** a published constant lies outside the certified
   interval — report the interval, the $d$, and the failing step.
   **Improve:** the certified interval is strictly tighter than published.
2. **Gate B (whole-chain enclosure table, hours–days).**
   Same treatment for Fertonani–Duman's four genie-aided bounds at
   $d\in\{0.25,0.5,0.64\}$ (including exact treatment of the infinite-series
   closed-form reduction, with a certified tail bound rather than a truncation),
   the Kanoria–Montanari 2- and 3-term small-$d$ coefficients in **exact
   rational** arithmetic, and Pinto–Ribeiro's $0.3578$. For Pinto–Ribeiro this
   doubles as the first CPU reproduction of a GPU-only result.
   **Pass/refutation criteria per row as in gate A.**
3. **Gate C (finite-length claims, days).**
   Independently recompute Morozov–Duman's finite-length achievability and
   converse bounds, and Tavakoli–Nguyen–Bose's q-ary sandwich at their own
   $(n,q)$ points, then extend to $(n,q)$ they did not run. A certified
   finite-blocklength improvement on their $d=1/2$ sandwich is the target.

## Pre-campaign requirements

- Owner first-hand reads of arXiv:0810.0785, 0912.5176, 1104.5546, 2604.05867,
  2504.20961, 2607.19559 (Rubinstein–Con is already owner-read). Extract per
  paper: the exact optimization solved, the alphabet/window truncation, the
  stated numerical precision, and the figure/table where each constant appears.
- Decide and commit the interval-arithmetic treatment of Blahut–Arimoto
  *before* running: BA is an iterative fixed-point method, so a naive interval
  iteration diverges. The pre-statement must fix the certified bounding
  strategy (e.g. dual/upper-bounding functional evaluated at a rounded
  primal iterate, with an explicit gap bound) and the truncation tail bound.
- One-page pre-statement per gate: precision, rounding mode, truncation order,
  tail bound, and numeric pass/fail thresholds.

## Disjointness

Classical Shannon theory. `../math/qec/` and all `../physics/` targets are
quantum; `../math/` has no information-theoretic target. No overlap.

## Layout

- `src/`, `campaigns/`, `scratch/` — created at first use.


## Current state (agent Delcap, 2026-08-30)

Status of each gate after this campaign (details in `campaigns/` rows and
`src/` headers). Campaign dirs are immutable; evidence labels per cs/README.md.

### Gate A (RC-2305.07156 crux rows) — PARTIAL, per protocol

* **A1/A2/A3/A4 (upper 0.3745(1-d) at d=0.68, 0.75, 0.85, 0.95): UNREACHED.**
  The published rows derive from Cbar_{28,k} for all k (U(28,0.68) at n=28,
  Table 1 entry 0.1199). Our certified schema (direct BA-dual certificates
  plus Arb-recurrence push of their Lemmas 2.3/2.4/2.5) certifies only to
  n<=17 base (direct) with recurrence to 28; the resulting certified
  U(28,0.68) = 0.1316873340 (+-1e-89) normalized /0.32 = 0.4115229, which is
  a valid certified UPPER bound but WEAKER than their 0.1199. Their numeric
  pipeline is simply out of reach on this workstation (memory/time wall:
  dense exact matrices for n>=18 with k>=11 exceeds 2^27 int64 cells; their
  sparse-schema near-k=n trick does not help the mid-k band that carries the
  d=0.68 binomial mass at n=28). No verdict offered either way on their
  constant. NOT presented as confirmation or refutation.
* **A5 (lower 0.1221(1-d)): REPRODUCED / COMPUTATIONAL-EVIDENCE.**
  [REPRODUCED] their exact public pipeline (BDC_Lower_Bounds repo, file
  hash list in campaigns/*prc_repro*), hyperparameters (0.19, 7.72/0.19,
  0.438), RZK table 1024x1024x128, step_limit=10000: rate
  0.12213624912193806 vs their notebook 0.12208469671934222 (paper rounds to
  0.1221). NOT a certificate: their k_probs clipping at 1e-300 makes the
  DM07 entropy term non-rigorous; the certified route would need their
  distribution input pinned under exact rational arithmetic (their PRC
  distribution is a floating-point BA output).

### Gate B (whole chain) — 12/12 FD10 Table II rows MACHINE-VERIFIED

* FD10 Table II (2-dp "upper bound" capacity entries f(L,R) for
  L<=7): all 12 published values ENCLOSED-BELOW-PUBLISH — our certified
  intervals sit strictly below every published figure, consistent with
  their stated round-up rule applied to a partially-converged BA DUAL (the
  dual approaches capacity from above). Offset is IRREGULAR (+0.0006 to
  +0.0127) and never a uniform +0.01: their published table is an upper
  bound catalog, not capacity estimates. Dissolved false positive (I
  initially suspected a +0.01 systematic shift; 12-pair measurement
  falsified that hypothesis). See `campaigns/.../direct_certificates.json`
  and `scratch/fd_pattern_report.txt`.
* **Falsified my own candidate reduction before building on it** (recorded
  in `scratch/`; log in `scratch/rc28_log.txt` lineage): a candidate
  S_n-equivariance "type-channel" capacity reduction was tested at
  n=3..7 and is FALSE (type-channel capacity is strictly below full
  capacity at every (n,k) tested except (3,1); e.g. (3,2): 1.4698 vs
  1.2539). Why: deletion preserves order and the output is a subsequence,
  so a permutation of the input DOES NOT induce a permutation of the
  output distribution (S_n-equivariance of the channel law fails). Not
  used anywhere in the pipeline; included here so the next person does not
  re-derive the same false premise.

### Gate C — NOT REACHED (runway)

Finite-length Morozov-Duman rows and the Tavakoli q-ary sandwich: not
started; needs the same certified-BA treatment at (n,q) grids. The base
pipeline (`src/delcap_cert.py`, `src/rc28.py`) is ready for it.

### The pipeline (deliverable independent of rows)

`src/delcap_cert.py` + `src/rc28.py`: certified-Arb Blahut-Arimoto.
mpmath 150-dps BA locating step -> exact rational snap (m/M on the simplex)
-> Arb outward-rounded evaluation of BOTH the primal I(p*) and the
Csiszar-Tusnady dual max_x KL(W(.|x) || D_{p*}) >= Cbar_{n,k}. Interval
[primal, dual] ALWAYS contains the true capacity regardless of BA
convergence; the BA only controls tightness. Empirically gap <= 1e-11 at
n<=4, ~3e-8 at n=8. Why this is sound (written in `src/delcap_cert.py`
header): (i) any full-support D' gives Cbar <= max_x KL (valid for ANY
BA iterate; the locating step's rounding errors only cost tightness, not
validity); (ii) every quantity evaluated is a log of an exact rational,
so the Arb ball arithmetic is one-shot outward-rounded evaluation of an
explicit formula — NO interval iteration anywhere, so nothing diverges.
Taint test recorded: ball widths grow monotonically 2.9e-180 -> 7.9e-90 ->
1.2e-35 -> 1.4e-17 as working precision drops 600 -> 60 bits (interval
path confirmed live at every stage, no float leak; float64 width checks
are degenerate and DO NOT detect this — see PROGRESS.md).

Full enclosure table (CSV) in `campaigns/enclosure_table.csv`.

---

## Gate C — DONE (agent DelcapGateC, 2026-08-30)

Finite-length claims for the binary/q-ary deletion channel. Gates A and B are
untouched; the ladder and pre-campaign requirements above are unchanged.

**Headline.** A certified finite-blocklength improvement on the
Tavakoli-Nguyen-Bose (arXiv:2607.19559) sandwich at the ticket's target
`d = 1/2`, and at their own `d` grid: **28 rows in which BOTH ends of the
published sandwich are strictly improved by a certified interval**, widths
`2.92e-10` to `6.73e-06` bits/symbol. Label `MACHINE-VERIFIED`.

**Frozen campaign** (scripts as run, raw outputs, deliverable CSVs, deviation
log, retraction record, tool versions, manifest, 39 checksums all verifying):
`campaigns/2026-08-30T13:07:27Z_7988b619-e89e-4e77-8de7-2c2a37a429d9_23445e8dbeeb/`
It indexes four source campaigns and regenerates its own anchor, second-route
and support-guard evidence at freeze time rather than transcribing them.

### First-hand reads (four extracted items each, recorded in `pre_statement.md`)

All three primaries were fetched and read from ar5iv full text by this agent;
none is `REPORTED`.

* **P1 arXiv:2504.20961** Morozov-Duman. Optimization: layer-oriented converse
  `M <= L(n,m,eps,Lambda) = (sum_{Lambda} tau_w)^n / ((sum_{Lambda} p_w)^n - eps)`,
  minimised over subsets `Lambda` of output-length classes (eqs 24-30).
  Truncation: complete `E(m,w)` integer tables for `m <= 23` (Table I), partial
  to `m <= 32` (Table II). Precision: `E(m,w)` exact integers; Table III code
  rates printed to 5 dp. Constants: Table I, eq (28), Table III, Fig. 2.
* **P2 arXiv:2604.05867** Pinto-Ribeiro. Optimization: GPU-parallelised
  Blahut-Arimoto on the exact deletion channel `C_{n,k}`, combined by the
  Fertonani-Duman convex combination (their Lemma 2). Truncation: all
  `k <= n <= 29`, plus `C_{31,k}` for `k <= 18`; no output-window truncation.
  Precision: BA tolerance `a = 0.005` (`a = 0.05` for `C_{31,k}`, `14<=k<=18`),
  reported values = BA rate + tolerance, an additive round-UP. Constants:
  their Tables 1, 2, 3.
* **P3 arXiv:2607.19559** Tavakoli-Nguyen-Bose. Optimization: the sandwich
  `LB1 = (1-d)log2 q - h2(d)`, refined
  `LB+ = (1-d)log2 q + H_Bin(n,1-d)/n - h2(d) + Delta_n(d)/n` (their Cor. 1),
  upper `UB = (1-d) log2 q` (their Thm 1), with
  `Delta_n(d) = sum_k w_k Phi_{k,n}` over exact integer pattern counts.
  Truncation: none — finite exact sums. Precision: printed to 3 dp in their
  Table I; `LB+` is an exact closed form. Constants: their Table I, 18 rows.

### Anchors reproduced before any new point

Two already-frozen gate-B rows, re-certified with **the Gate C certificate code
path itself** (`gate_c_dhalf.cert_primal` / `cert_dual`), not merely with the
older pipeline — so the check exercises the layer that produced every Gate C
number. The channel-builder layer is anchored separately and by different
evidence: the Tavakoli Examples 1-3 and Thm 5 self-test (exact pattern
histogram `{1:4, 2:2}`, `Phi_{1,2}(q)=1/q`, `Delta_2(1/2)=1/4`).

| frozen row | frozen printed interval | Gate C code path | width | contained |
|---|---|---|---|---|
| `f(3,2)` | `[1.4697819938, 1.4697820261]` | `[1.4697819937559, 1.4697819937572]` | `1.3627e-12` | yes |
| `f(5,3)` | `[1.8715442105, 1.8715442453]` | `[1.8715442104903, 1.8715442104926]` | `2.2151e-12` | yes |

Both certified intervals are contained in the frozen ones (the frozen endpoints
are printed to 10 dp, and both of our endpoints agree with them at that
precision). The `f(3,2)` upper endpoint additionally reproduces the
*independent* frozen `rc28.py` direct dual certificate `1.469781993757231` with
`|ours - rc28| = 0.000e+00` at float resolution — a second route to the same
endpoint from different code, different snap and a different locator. Raw
output: `anchor_check.txt` in the freeze campaign, regenerated by the freeze
script itself rather than transcribed.

### Second route for the `LB+` column (overdetermination, not trust)

The `their LB+` column below carries the conclusion, so it is computed twice by
machinery that shares nothing. Route 1 is their closed form (their Cor. 1):
`h2(d)`, `H_Bin(n,1-d)` and `Delta_n(d) = sum_k w_k Phi_{k,n}` from exact
integer pattern counts. Route 2 is the certified exact mutual information at
**uniform** input, evaluated from the full channel transition matrix by
`cert_primal` — a different formula over different data structures. Their
Cor. 1 asserts these are the same number; that assertion is now
`MACHINE-VERIFIED` at seven points spanning `q` in `{2,3,4}`, `n` in `{3,5}`,
`d` in `{1/2, 1/20, 1/10, 1/5}`, agreeing to the last certified digit
(`|diff| = 0.00e+00` in every case, both routes' intervals coinciding). This
simultaneously validates the channel builder, `cert_primal`, and the sandwich
implementation, because a defect in any one of them would break the agreement.

### Target result: certified improvement on the d = 1/2 sandwich

Per-symbol. `LB+` is their tightened lower bound (their Cor. 1, the exact
uniform-input rate), `UB = (1-d) log2 q` their upper. Truncation order: **none**
— the full exact `q^n x sum_k q^k` deletion channel is enumerated, so the tail
bound is **exactly 0**. Arb working precision 400 bits; snap denominators
`2^30` (input) and `2^30`/`2^40` (output).

| q | n | their LB+ | their UB | certified interval | width | lower gain | upper gain |
|---|---|---|---|---|---|---|---|
| 2 | 2 | 0.37500000 | 0.50000000 | `[0.415241012, 0.415241013]` | `6.72e-10` | +0.0402 | -0.0848 |
| 2 | 3 | 0.30698620 | 0.50000000 | `[0.365269059, 0.365269061]` | `2.22e-09` | +0.0583 | -0.1347 |
| 2 | 4 | 0.26415061 | 0.50000000 | `[0.332341253, 0.332341258]` | `5.28e-09` | +0.0682 | -0.1677 |
| 2 | 5 | 0.23463455 | 0.50000000 | `[0.308784180, 0.308784202]` | `2.19e-08` | +0.0741 | -0.1912 |
| 2 | 6 | 0.21300111 | 0.50000000 | `[0.290941571, 0.290941731]` | `1.59e-07` | +0.0779 | -0.2091 |
| 2 | 7 | 0.19641139 | 0.50000000 | `[0.276852922, 0.276853161]` | `2.39e-07` | +0.0804 | -0.2231 |
| 2 | 8 | 0.18324295 | 0.50000000 | `[0.265375771, 0.265382499]` | `6.73e-06` | +0.0821 | -0.2346 |
| 3 | 2 | 0.62581458 | 0.79248125 | `[0.667481250, 0.667481251]` | `1.01e-09` | +0.0417 | -0.1250 |
| 3 | 3 | 0.53285995 | 0.79248125 | `[0.594406225, 0.594406261]` | `3.60e-08` | +0.0615 | -0.1981 |
| 3 | 4 | 0.47325262 | 0.79248125 | `[0.546348770, 0.546348925]` | `1.54e-07` | +0.0731 | -0.2461 |
| 3 | 5 | 0.43155955 | 0.79248125 | `[0.512042240, 0.512043064]` | `8.24e-07` | +0.0805 | -0.2804 |
| 4 | 3 | 0.70664941 | 1.00000000 | `[0.764131514, 0.764131527]` | `1.31e-08` | +0.0575 | -0.2359 |
| 4 | 4 | 0.63826289 | 1.00000000 | `[0.707114444, 0.707114895]` | `4.51e-07` | +0.0689 | -0.2929 |

Verdict on every row: `CERT_LOWER_BEATS_LBplus + CERT_UPPER_BEATS_UB`, decided
by Arb-endpoint inequalities (certified lower endpoint vs the exact-Arb `LB+`
upper endpoint), with margins 0.04-0.08 (lower) and 0.08-0.29 (upper)
bits/symbol — four to six orders of magnitude above any interval width.

### Certified recomputation at the paper's own (q,n,d) points

| q | n | d | their LB+ | their UB | certified interval | width | lower gain | upper gain |
|---|---|---|---|---|---|---|---|---|
| 2 | 3 | 1/20 | 0.90976488 | 0.95000000 | `[0.910574224, 0.910574225]` | `2.92e-10` | +0.0008 | -0.0394 |
| 2 | 3 | 1/10 | 0.82451503 | 0.90000000 | `[0.827679159, 0.827679161]` | `2.14e-09` | +0.0032 | -0.0723 |
| 2 | 3 | 1/5 | 0.66847117 | 0.80000000 | `[0.680473561, 0.680473562]` | `1.24e-09` | +0.0120 | -0.1195 |
| 2 | 5 | 1/20 | 0.88619813 | 0.95000000 | `[0.887557899, 0.887557903]` | `3.36e-09` | +0.0014 | -0.0624 |
| 2 | 5 | 1/10 | 0.78229360 | 0.90000000 | `[0.787495955, 0.787495959]` | `4.10e-09` | +0.0052 | -0.1125 |
| 2 | 5 | 1/5 | 0.60155006 | 0.80000000 | `[0.620298060, 0.620298068]` | `7.85e-09` | +0.0187 | -0.1797 |
| 2 | 10 | 1/20 | 0.85098805 | 0.95000000 | `[0.852997769, 0.852997906]` | `1.37e-07` | +0.0020 | -0.0970 |
| 2 | 10 | 1/10 | 0.72247657 | 0.90000000 | `[0.730086455, 0.730086624]` | `1.69e-07` | +0.0076 | -0.1699 |
| 2 | 10 | 1/5 | 0.51573107 | 0.80000000 | `[0.542096537, 0.542096780]` | `2.42e-07` | +0.0264 | -0.2579 |
| 3 | 3 | 1/20 | 1.45321966 | 1.50571438 | `[1.453878938, 1.453878943]` | `5.07e-09` | +0.0007 | -0.0518 |
| 3 | 3 | 1/10 | 1.32766925 | 1.42646625 | `[1.330307006, 1.330307014]` | `7.85e-09` | +0.0026 | -0.0962 |
| 3 | 3 | 1/5 | 1.09470126 | 1.26797000 | `[1.105231693, 1.105231699]` | `6.27e-09` | +0.0105 | -0.1627 |
| 3 | 5 | 1/20 | 1.42508255 | 1.50571438 | `[1.426057398, 1.426057434]` | `3.58e-08` | +0.0010 | -0.0797 |
| 3 | 5 | 1/10 | 1.27653927 | 1.42646625 | `[1.280432698, 1.280432727]` | `2.82e-08` | +0.0039 | -0.1460 |
| 3 | 5 | 1/5 | 1.01112755 | 1.26797000 | `[1.026499717, 1.026499762]` | `4.50e-08` | +0.0154 | -0.2415 |

Agreement is explicit: every certified interval lies strictly inside the
published sandwich, so all 15 rows **agree** with P3 and none excludes a
published value. Their `LB+` is reproduced as a lower endpoint to `1e-8` at
`n=3` before our BA improvement is applied.

**Swept set, exactly (rule 7).** `d = 1/2`: `(q,n)` in `(2,2..8)`, `(3,2..5)`,
`(4,3)`, `(4,4)`. Paper grid: `(q,n)` in `(2,3)`, `(2,5)`, `(2,10)`, `(3,3)`,
`(3,5)` at each `d` in `{1/20, 1/10, 1/5}`. **NOT swept:** `(3,10)` at any `d`
(exact channel is `3^10 x 88573`, beyond this workstation) and the
pre-registered `d=1/2` extensions `(3,7)`, `(3,8)`. This is a finite-`n`
statement about `C_{q,n}(d)`; it says nothing about the `n -> infinity`
capacity.

### Morozov-Duman Table III: 36 rows recomputed and certified

All 36 LO-CVB rows (`m` in `{5,22,23}`, `n` in `{1,...,1024, infinity}`,
`delta = eps = 1/5`) independently recomputed from the integer `E(m,w)` tables
and certified: interval widths `0.0` at every printed digit, maximum ball
radius `2.84e-120` over all 36 rows. Truncation order: none (finite closed
form); tail bound exactly 0. Label `MACHINE-VERIFIED` for our values.

The published values agree with ours to within **`7.1e-7` absolute in every
row** under the paper's evident round-up printing convention, which is the
correct conservative choice for a converse bound. Precisely:
`gap = printed - certified` lies in `[7.10008e-07, 1.03994e-05]` with spread
`9.689365e-06`, strictly **less than one 5-dp ulp**; solving
`printed = ceil_5dp(certified + Delta)` over all 36 rows simultaneously gives
`Delta` in `(3.990e-07, 7.10008e-07]`. That band is **determined, not fitted**
— under ceiling printing `gap` lies in `[Delta, Delta+ulp)`, so the minimum gap
caps `Delta` above and the maximum gap bounds it below — and its non-emptiness
is *equivalent* to spread < 1 ulp. Under round-half-up the band is
`[5.399e-06, 5.710e-06]`.

This is a **precision-level effect, not a formula difference**: a differing
formula term would be `n`- or `m`-dependent and would break a single-`Delta`
fit spanning `n = 1` to `n = infinity` and `m = 5, 22, 23`. The fit holding
across all 36 rows is positive evidence against a formula difference. The
`~1e-6` relative size is that of a 6-to-7-significant-digit intermediate.

**AMENDED 2026-08-30, and the amendment is recorded rather than substituted.**
This section first characterised the offset as a *resolved* round-up printing
convention. That characterisation is **withdrawn**. The reason is the TNB
Table I test below, which I ran precisely to check whether conservative
round-up printing is a convention in this literature: it is not. Two of the
three primaries print conservative round-ups (Morozov-Duman by measurement
here, Pinto-Ribeiro by their own stated `BA rate + tolerance`), but
Tavakoli-Nguyen-Bose round to nearest, inconsistently, with no conservative
bias. So the offset does **not** become typical of the field, and it is now
recorded as an **unexplained residual** with its direction and magnitude
stated: printed exceeds certified in all 36 rows, by `7.10008e-07` to
`1.03994e-05`. A single constant `Delta` reproduces every printed value, which
bounds the effect sharply, but nothing here explains its cause. The earlier
`~1e-6`-relative mechanism note remains `INFERENCE` and is **not** strengthened
by the TNB test.

These are the first certified enclosures of those rows. **No published claim is
contradicted; their converse stands.** Escalated to Main before anything was
written, and cleared by Main after independent re-derivation of the band.

`E(m,w)` provenance: `m=5` all `w` recomputed from scratch by brute force,
reproducing their eq (28) values `E(5,2)=32, E(5,3)=52, E(5,4)=54` exactly
(`MACHINE-VERIFIED`). For `m=20..23` the columns `w` in `{0,1,2,m}` are
re-derived independently and agree with their Table I — `w=0`: 1; `w=1`: `2m`;
`w=m`: `2^m`; `w=2`: closed form `2*C(m,2) + 2*floor(m^2/4)`, proved in
`src/gate_c_md.py` and brute-force validated at `m=4..12`, reproducing their
`580, 640, 704, 770` (`MACHINE-VERIFIED` for those cells). Remaining cells:
`CITED-DEPENDENCY`. Aggregate check: the `n->infinity` row is a weighted sum
over the whole `E(m,.)` column and reproduces `0.80272`, `0.73569`, `0.73414`.

### Tavakoli-Nguyen-Bose Table I: all 18 printed entries transcribed and tested

Their Table I was read first-hand and all 18 printed entries transcribed, then
compared against the exact-Arb recomputation of their own Cor. 1 / Thm. 1
closed forms. This was run to answer a specific question: two of the three
primaries in this target print conservative round-ups, so **is conservative
round-up printing a convention in this literature?**

**Answer: no.** Printed minus computed, with signs:

| column | rows | distinct | printed above | below | exact | diff range |
|---|---|---|---|---|---|---|
| `LB1` (`n`-independent) | 18 | 6 | 3 | 15 | 0 | `[-4.71e-4, +3.97e-4]` |
| `LB+` | 15 | 15 | 7 | 8 | 0 | `[-5.39e-4, +4.85e-4]` |
| `UB` (`n`-independent) | 18 | 6 | 6 | 3 | 9 | `[-4.66e-4, +2.86e-4]` |

For the `LB+` column **no printing convention is feasible**: solving
`printed = round(computed + Delta)` at 3 dp over all rows simultaneously gives
an EMPTY band under round-half-up, under ceiling and under truncation. The
reason is structural and is the exact mirror of the Morozov-Duman case, under
the same criterion — a single-offset fit exists **iff** the spread of
printed-minus-computed is below one ulp:

* MD Table III at 5 dp: spread `9.6890e-06` **<** `1.0e-05` ulp → `Delta` exists
* TNB Table I `LB+` at 3 dp: spread `1.0240e-03` **>** `1.0e-03` ulp → none exists

What their digits are instead: signs mixed nearly evenly and every `|diff|` at
most `5.39e-4`, i.e. about half an ulp — correct rounding, applied
inconsistently row to row. Two rows in one table under two different rules:
computed `0.824515` printed `0.825` is reachable ONLY by round-to-nearest
(truncation gives `0.824`), while computed `1.276539` printed `1.276` is
reachable ONLY by truncation (round-to-nearest gives `1.277`).

Verdict on their formula columns: **agreement at their printed precision** —
`LB1` and `UB` reproduced on 18 of 18 rows, `LB+` on 15 of 15, every deviation
within half an ulp. The three `q=3, n=10` `LB+` entries are transcribed and
labelled not-recomputed (`Delta_n` needs the pattern counts at `3^10`, the
pre-registered drop).

### Certified lower bounds above the published BA column

Their `C_{(q,n)}` column is their Blahut-Arimoto value for the finite-block
capacity, printed to 3 dp. Our certified lower endpoint is the exact mutual
information of an explicit rational input distribution in outward-rounded Arb
at 400 bits, so this comparison rests on the **primal certificate alone** and
does not involve the dual.

**In 9 of the 15 compared rows we now hold a tighter and higher certified lower
bound on the finite-block capacity than the published BA column**, by more than
one 3-dp ulp:

| q | n | d | their `C_{(q,n)}` | our certified interval | our lower bound exceeds theirs by |
|---|---|---|---|---|---|
| 2 | 3 | 1/5 | 0.676 | `[0.680473561, 0.680473562]` | `+0.004474` |
| 2 | 5 | 1/10 | 0.786 | `[0.787495955, 0.787495959]` | `+0.001496` |
| 2 | 5 | 1/5 | 0.613 | `[0.620298060, 0.620298068]` | `+0.007298` |
| 2 | 10 | 1/10 | 0.728 | `[0.730086455, 0.730086624]` | `+0.002086` |
| 2 | 10 | 1/5 | 0.531 | `[0.542096537, 0.542096780]` | `+0.011097` |
| 3 | 3 | 1/10 | 1.329 | `[1.330307006, 1.330307014]` | `+0.001307` |
| 3 | 3 | 1/5 | 1.101 | `[1.105231693, 1.105231699]` | `+0.004232` |
| 3 | 5 | 1/10 | 1.279 | `[1.280432698, 1.280432727]` | `+0.001433` |
| 3 | 5 | 1/5 | 1.020 | `[1.026499717, 1.026499762]` | `+0.006500` |

Five further rows exceed by less than one ulp (`5.7e-05` to `9.98e-04`), and one
row (`q=3, n=3, d=1/20`, printed `1.454`) sits above our interval consistent
with round-to-nearest.

**Their sandwich is confirmed, not challenged.** Their stated ordering is
`LB1 <= LB2 <= LB+ <= C_{q,n} <= UB`, and our certified interval lies strictly
inside their own `[LB+, UB]` in all 15 compared rows. Their Section V claims
that the BA values exceed `LB+`, which their printed numbers do.

Reading of their column, labelled `INFERENCE` and resting on structure rather
than on any single row: it is an **unconverged BA estimate** rather than a claim
of the capacity to printed precision. The excess is monotone increasing in `d`
at fixed `n` and monotone increasing in `n` at fixed `d` (`q=2, d=1/5`:
`+4.474e-03`, `+7.298e-03`, `+1.1097e-02` at `n = 3, 5, 10`). Blahut-Arimoto's
primal rate climbs toward capacity **from below**, so an unconverged run
understates, and a certified lower bound exceeding it is what that predicts.
Our intervals are `MACHINE-VERIFIED`; the direction and magnitudes are exact.

### What printed tables in this literature are worth, measured

Across the three primaries this target consumes, printed numeric tables differ
from certified recomputation for **four distinct reasons**, and not one of them
touches a theorem:

1. **Morozov-Duman Table III** — conservative round-up: a single
   `Delta <= 7.10008e-07` plus 5-dp ceiling reproduces all 36 rows; cause of the
   offset itself unexplained.
2. **Pinto-Ribeiro** — author-stated additive round-up, `BA rate + tolerance`,
   which makes their printed numbers valid upper bounds by construction.
3. **TNB Table I formula columns** — inconsistent correct-rounding; no single
   `Delta` exists; every deviation within half an ulp.
4. **TNB `C_{(q,n)}` column** — an unconverged BA estimate, excess monotone in
   both `d` and `n`.

The measured observation: in finite-length deletion-channel work, printed tables
are not reliable to their printed precision, for several independent reasons.
Certified recomputation is therefore not redundant with reading the paper.

### New certified rows with no published analogue

* **30 LO-CVB extension rows** at `delta` in `{1/20, 1/2, 4/5}`, `m` in
  `{22,23}`, `n` in `{1,4,64,1024,infinity}`, `eps=1/5` — Morozov-Duman show
  these `delta` only as Fig. 2 curves, so there is no printed value to compare.
  Widths `0.0` at every printed digit. `MACHINE-VERIFIED`.
* **30 Pinto-Ribeiro `C_{n,k}` extension rows**, `n = 6,7,8,9` and all
  `1 <= k <= n`. Pinto-Ribeiro publish `C_{n,k}` only at `n in {29,31}` (GPU),
  so these are new certified values, not a reproduction, and their published
  rows are out of certified reach here. Because there is no published value to
  check them against, each row was produced by **two routes** that share only
  the exact integer channel matrix: route 1 float64 locator -> snap `2^30`/`2^40`
  -> `gate_c_dhalf` certificates at 400 bits; route 2 the frozen `delcap_cert`
  mpmath-150dps locator -> snap `2^24` -> `certify_ab` at 300 bits. Both
  intervals must contain the same `C_{n,k}`, so the test is overlap — a
  disjoint pair would prove a certificate wrong. **24 of the 30 rows carry both
  routes and all 24 overlap; 0 disjoint.** The remaining 6 (`n=9`, `k=3..8`)
  are route-1 only and labelled `FASTPATH_ONLY_no_mpmath_row`, because an
  mpmath row at `n=9` costs 10-80 min against ~1 s for route 1. Widths
  `0.0e+00` to `4.9e-04`. Full table with both routes' endpoints:
  `TABLE_pr_cnk_rows.csv` and `pr_extra_mpmath_crosscheck.txt`.
* **Greedy achievability (their Alg. 1) at `m=5`, `delta=1/5`**: exact rational
  frame-error rate for code sizes `M = 1..32` (e.g. `M=4` FER `0.023920`,
  `M=8` FER `0.138200`), computed in exact integer arithmetic.

### Solver / tolerance audit, and one live invalid-certificate branch caught

Prompted by a repo-wide warning: no LP/QP/MILP anywhere in the Gate C certified
path (grep for `scipy`/`linprog`/`highs`/`optimize` across all four gate-C
scripts returns 0). Float appears only in the BA locating step and in ranking
`Lambda` candidates; both feed an exact-rational snap.

The audit found a real hazard and confirmed the code is safe against it. The
dual certificate is `max_x KL(W(.|x) || D')`; if `D'(y) = 0` where
`W(y|x) > 0`, that KL is `+infinity`, and **silently skipping such a term would
produce a too-small "upper bound" — an invalid certificate that looks clean**.
Verified at `q=2, n=10, d=1/20`: the BA output marginal's minimum mass is
`9.765625e-14` while the snap resolution `2^-40 = 9.094947e-13` is `9.3x`
coarser, so exactly 1 of 2047 entries snaps to zero **and that output is
reachable**. `cert_dual` returns `+inf` (no claim) on that branch rather than
skipping; full support is additionally guaranteed by a `+1` bump and by an
always-full-support uniform `D'` candidate, taking the min (`8.529977690` vs
`9.672330910`). Same failure shape as a tolerance coarser than the quantity
being resolved, caught before it bit.

### Retractions (rule 5, kept inline)

* **R1 RETRACTED — "the LO-CVB rows disagree with the paper."** My original
  escalation said 19 of 36 rows "disagree", cause "last-digit rounding in a
  printed table". Falsified twice over: (a) Main showed by exact decimal
  arithmetic that my unshifted certified values round to `0.66390`/`0.80278`
  under both round-half-up and ceiling while the printed values are
  `0.66391`/`0.80279`, differences `1.040e-5`/`1.037e-5` exceeding the `1.0e-5`
  ulp, so a rounding artifact of my own value cannot explain it; (b) I then
  falsified my own 19/17 partition — sorted by gap it splits exactly at my
  self-chosen `+-5e-6` window, carries zero structural information, and the
  classes are interleaved on every property (`m` on both sides, full-support
  `Lambda` 11/19 vs 6/17, `n` from 1 to infinity on both sides). Gap
  correlations are noise at `n=36`: `n` +0.40, `m*n` +0.28, `is_full_set`
  +0.20, `m` -0.08, `|Lambda|` -0.03. A third hypothesis — a systematic
  additive finite-`n` term — was also falsified: the gaps are a continuum, 35
  distinct values among 36 rows at `1e-8` resolution. **Correction I
  volunteered against my own report: the direction split is 36-0, not 19-0.**
  Replaced by the measured statement in that section.
* **A2 AMENDED — the Morozov-Duman offset is no longer called a resolved
  printing convention.** The amendment sits inline in the Morozov-Duman section
  above, where the original characterisation was made. Cause: the TNB Table I
  test, which I ran to check whether conservative round-up printing is a
  convention in this literature, showed it is not — so the offset does not
  become typical and is now recorded as an unexplained residual with its
  direction and magnitude stated. Discharging a requirement destroyed the
  hypothesis the discharge was meant to support; that is recorded, not hidden.
* No other Gate C claim has been retracted.

### Deviations from the committed grid

`pre_statement.md` was committed before the first run and is **not** edited
retroactively; every deviation is logged in `DEVIATIONS.md` in the freeze
campaign. The material one: the pre-statement fixes an mpmath-150dps locating
step, and these rows used a float64 locator instead (the mpmath locator cost
~290 s/row at `q=3,n=5` and collapsed onto a near-degenerate `p` under heavy
deletion). This cannot affect validity — the certified lower bound is a true
mutual information of an exact rational `p*` and the certified upper bound is a
true `max_x KL` against an exact rational full-support `D'`, for **any**
candidates whatsoever — it affects only width, and the achieved widths are
tighter than the mpmath path's. Also dropped: `(3,10)`, `(3,7)`, `(3,8)`,
`n in {10,11,12}` for Pinto-Ribeiro, and `m in {10,15}` for Morozov-Duman
(`m=22` substituted, being a column their Table III actually prints).

### What remains open

* `(q,n) = (3,10)` and the `q=3` ladder above `n=5` need either a sparse
  channel representation or a GPU; the dense exact-integer matrix is the wall.
* The Morozov-Duman Table III offset is an **unexplained residual**: printed
  exceeds certified in all 36 rows by `7.10008e-07` to `1.03994e-05`, a single
  constant `Delta` reproduces every printed value, and the field-convention
  explanation is falsified by the TNB test. Its cause is open. Settling it needs
  their intermediate `tau`/`p` values or their code at higher precision, neither
  of which is available here, and it will not be guessed at.
* The `Lambda` search at `m=22,23` certifies only the top 24 float-ranked
  candidates per row, so the reported LO-CVB is not a *proven* global optimum
  over `Lambda` (it remains a valid bound either way).
* Pinto-Ribeiro's published `n in {29,31}` rows remain out of certified reach
  here, so their headline `0.3578(1-d)` is untouched by this gate.
* Nothing here bears on the `n -> infinity` capacity `C(BDC_d)`, which stays
  open; every Gate C row is a finite-`n` theorem.



## Current state — Open Item 1 resolved: sparse orbit certificates for
## (3,10) and the q=3 ladder (agent DelcapNextGate, 2026-08-30)

**Outcome first.** The first listed open item — "`(q,n) = (3,10)` and the
`q=3` ladder above `n=5` … the dense exact-integer matrix is the wall" — is
**resolved by a sparse group-orbit certificate representation**, on CPU, with
no GPU: **20 pre-registered rows certified, and in ALL 20 both ends of the
published Tavakoli-Nguyen-Bose sandwich are strictly improved by certified
intervals**, label `MACHINE-VERIFIED`, widths `9.9e-08` to `1.32e-04`
bits/symbol. The wall row `(3,10,d=1/2)` is certified per-symbol
`[0.42369976, 0.42371229]`, width `1.25e-05`, inside their sandwich
`[0.32826401, 0.79248125]` with margin 0.095 on the lower end and 0.369 on
the upper — three to five orders of magnitude above the interval width.
Frozen campaign:
`campaigns/2026-08-30T17:50:25Z_ddd57c90-d088-47c1-acda-dcb54bcf3555/`
(pre-statement committed before the first computation; `manifest.md`;
scripts byte-frozen; 9-file sha256 set all verifying; two mid-run bugs and
their fixes recorded inline in the manifest).

### The reduction, and how it avoids the falsified gate-B route

The group is `G = S_q x C_2` acting on **symbol VALUES** (alphabet relabeling
`sigma` plus word reversal `rho`), applied simultaneously to input and output
words. Deletion acts on **positions**; value edits commute with position
deletion, so `A(sigma x, sigma y) = A(x,y)` and `A(rho x, rho y) = A(x,y)`
with `A(x,y)` = number of subsequences of `x` equal to `y`. This is NOT the
`S_n` input-POSITION/type reduction falsified in gate B (deletion preserves
order, so position permutations are not a symmetry); Main independently
re-derived both equivariances in exact arithmetic at `(3,2,1/3), (4,2,1/5),
(3,3,2/7)` and confirmed the position-permutation control BREAKS — the
falsified reduction and this one are separated by measurement on both sides.
The pipeline additionally anchor-tests the equivariance identity pointwise
at runtime (exhaustive at `q=2,n=3`; KL-invariance with a non-invariant
control that the test can detect, max deviation `7.5e-01`), plus orbit-size
closure, exact p-marginal agreement with the dense word enumeration
(`|diff| = 0.00e+00`), Phi-route agreement (`<= 1e-12` at 5 points), and
primal exactness vs the dense word-level evaluation (`<= 1e-12`, uniform +
skew inputs). Dense-vs-orbit cross-validation at `(3,6),(3,7) d=1/2`:
**OVERLAP OK**, orbit lower endpoint matching the dense lower endpoint to 14
digits (`crosscheck.json`).

Certificates (all exact rationals end to end, Arb 400-bit outward rounding):
primal = exact mutual information of one explicit rational input `p*`
(valid for any `p*`; `D'` derived as the true marginal of `p*` from merged
orbit dictionaries); dual = `max_x KL(W(.|x) || D')` over input g-orbit
representatives for exact rational `G`-invariant `D'` candidates — validity
of the Csiszar-Tusnady upper bound needs only full support, unchanged from
the frozen pipeline; the orbit identities shrink evaluation cost only. The
frozen `cert_dual = +infinity` no-claim branch is preserved and re-verified
VERBATIM at `q=2, n=10, d=1/20` (raw per-word snap at `2^-40` zeros exactly 1
reachable word, mass `9.766e-14 < 9.095e-13`; `cert_dual` returns `+inf`).
The float64 BA locator runs in the invariant simplex with per-rep sparse
rows (~750k nonzeros at `(3,10)` vs `5.2e9` dense cells) and feeds only the
snap. The dense wall is not climbed; it is bypassed.

### The 20 certified rows (finite-`n` theorems about `C_{3,n}(d)`)

Per-symbol. `their LB+` / `their UB` from arXiv:2607.19559 Cor. 1 / Thm. 1
recomputed exactly in Arb (their Table I prints `(3,10)` to 3 dp; our
exact-Arb LB+ agrees with those prints at `(3,10)` — `1.386->1.38574`,
`1.208->1.20833`, `0.908->0.90848` — inside their 3-dp rounding).

| n | d | their LB+ | their UB | certified interval | width |
|---|---|---|---|---|---|
| 6 | 1/2 | 0.40060333 | 0.79248125 | `[0.48612691, 0.48613104]` | `4.13e-06` |
| 6 | 1/5 | 0.98243193 | 1.26797000 | `[0.99928060, 0.99928076]` | `1.62e-07` |
| 6 | 1/10 | 1.25818150 | 1.42646625 | `[1.26245154, 1.26245165]` | `1.07e-07` |
| 6 | 1/20 | 1.41473940 | 1.50571438 | `[1.41580723, 1.41580733]` | `9.86e-08` |
| 7 | 1/2 | 0.37659804 | 0.79248125 | `[0.46572597, 0.46575329]` | `2.73e-05` |
| 7 | 1/5 | 0.95898267 | 1.26797000 | `[0.97695052, 0.97695090]` | `3.82e-07` |
| 7 | 1/10 | 1.24280756 | 1.42646625 | `[1.24735973, 1.24735999]` | `2.65e-07` |
| 7 | 1/20 | 1.40595463 | 1.50571438 | `[1.40709141, 1.40709165]` | `2.39e-07` |
| 8 | 1/2 | 0.35735985 | 0.79248125 | `[0.44916027, 0.44929226]` | `1.32e-04` |
| 8 | 1/5 | 0.93940600 | 1.26797000 | `[0.95824532, 0.95824676]` | `1.43e-06` |
| 8 | 1/10 | 1.22968245 | 1.42646625 | `[1.23445209, 1.23445306]` | `9.66e-07` |
| 8 | 1/20 | 1.39835199 | 1.50571438 | `[1.39954151, 1.39954245]` | `9.44e-07` |
| 9 | 1/2 | 0.34154094 | 0.79248125 | `[0.43538141, 0.43538579]` | `4.38e-06` |
| 9 | 1/5 | 0.92278559 | 1.26797000 | `[0.94231948, 0.94232456]` | `5.07e-06` |
| 9 | 1/10 | 1.21830837 | 1.42646625 | `[1.22324996, 1.22325236]` | `2.40e-06` |
| 9 | 1/20 | 1.39167586 | 1.50571438 | `[1.39290682, 1.39290944]` | `2.61e-06` |
| 10 | 1/2 | 0.32826401 | 0.79248125 | `[0.42369976, 0.42371229]` | `1.25e-05` |
| 10 | 1/5 | 0.90848105 | 1.26797000 | `[0.92857975, 0.92859610]` | `1.63e-05` |
| 10 | 1/10 | 1.20833242 | 1.42646625 | `[1.21341296, 1.21342317]` | `1.02e-05` |
| 10 | 1/20 | 1.38574454 | 1.50571438 | `[1.38700885, 1.38701617]` | `7.32e-06` |

Verdict on every row: `CERT_LOWER_BEATS_LBplus + CERT_UPPER_BEATS_UB`,
decided by Arb-endpoint inequalities at 400 bits — the certified lower
endpoint exceeds their exact `LB+` upper endpoint (by `+6.6e-03` at
`(6,1/20)` up to `+9.5e-02` at `(10,1/2)`), and the certified upper endpoint
sits below their exact `UB` lower endpoint (by `0.21` at `(10,1/2)` to
`0.049` at `(10,1/20)`). All margins are three to six orders of magnitude
above the interval widths. **No published value is excluded by any interval:**
every certified interval lies strictly inside their own `[LB+, UB]`, so this
confirms their sandwich theorem while tightening both ends at their own
finite-`n` points — the same verdict shape as the gate-C `d=1/2` campaign.
Their BA column at `(3,10)` prints `1.387 / 1.212 / 0.920`; our certified
lower endpoints `1.38700885 / 1.21341296 / 0.92857975` exceed the first two
by under one 3-dp ulp and the third by `+8.6e-03` — the same
unconverged-BA-estimate reading recorded at gate C (`INFERENCE`; the
intervals are `MACHINE-VERIFIED`).

### Rule-7 sweep sentence (mandatory)

**Swept set exactly:** `q = 3`; `(q,n,d)` = `(3,10)` at `d` in
`{1/2, 1/5, 1/10, 1/20}` and `(3,n)` at `n` in `{6,7,8,9}` at the same four
`d` — 20 rows certified (plus the two cross-validation rows at `(3,6),(3,7)
d=1/2` whose dense-path twins come from the frozen pipeline); Arb 400-bit
outward rounding, exact `fmpq` rationals, snap denominators input `2^30` /
output per-word `2^30`, float64 locator feeding only the snap. **NOT
swept:** `q >= 4` at any `n > 5`; any `n > 10`; Morozov-Duman `m > 23`; the
Lambda-search optimality gap at `m` in `{22,23}`; Pinto-Ribeiro `n` in
`{29,31}`; and the `n -> infinity` capacity `C(BDC_d)` — every row above is
a finite-`n` theorem about `C_{3,n}(d)` and is asserted as nothing more
(rules 7/9/14); the asymptotic side of this target remains open exactly as
recorded above.

### What of the open item remains open

The sparse representation removes the dense-matrix wall for `q=3` through
`n=10` and, with the same machinery, makes `n=11,12` and `q=4` ladders
mechanically reachable (not pre-registered in this campaign's box; a
separate pre-registration would be required before running them — rule 16).
`C_{3,10}(d)` now has certified two-sided brackets at four `d` values.
Escalation to Main: none triggered — no interval excluded any published
value.

## Owner correction — 2026-08-31: prior $q=3$ representative certificates suspended

The “20 certified rows” claim immediately above is **RETRACTED pending
orbit-safe replay**. A mandatory full-alphabet audit found that the frozen
largest-remainder snap need not preserve input-orbit invariance, although the
subsequent dual maximization compressed to one representative per orbit.

The first audited row, $(q,n,d)=(3,6,1/2)$, passed: exact orbit equality,
$3{,}279$ generator checks, all $729$ input divergences, zero
full-versus-representative overlap failures, and identical representative/full
maxima. The second row, $(3,6,1/5)$, failed **before** any new bound was
accepted: the frozen `ba_word` integer vector at denominator $2^{30}$ split
flat output orbit 51 (size 12) into counts $578{,}028$ and $578{,}029$.
Therefore that candidate's output law is not constant on output orbits, and a
representative-only maximum is not a valid upper certificate for it.

Audit status: **1/20 re-certified, 1/20 failed, 18/20 not run** under the
mandatory stop rule. Until a new campaign constructs probabilities in exact
input-orbit coordinates and replays the full alphabet, only the
$(3,6,1/2)$ row retains certification from this table. The other 19 displayed
intervals remain preserved as historical outputs but are **UNVERIFIED /
SUSPENDED**, not machine-certified results. This is a repository-certificate
defect; it neither refutes nor changes any Tavakoli–Nguyen–Bose theorem or
published bound.

The proving counterexample is independently visible already at
$(q,n,d)=(3,3,1/2)$: a split snap makes the inadmissible representative-only
dual strictly smaller than the true full-alphabet maximum. The replacement
must snap orbit masses first, distribute each exact mass uniformly within its
orbit, assert exact output-orbit invariance, and compare the compressed
maximum against every full-alphabet input before any row can be promoted.

**Correction to the correction:** the demonstrated split is in the
$1{,}093$-entry **output reference distribution $D$**, not in the input
distribution $p$. The dual may use any positive normalized $D$; its safe
replacement is to snap total **output-orbit** masses and set
$D_y=D(O)/|O|$ exactly (or to induce $D$ exactly from an invariant input).
Input-orbit snapping is required only when the primal itself is compressed by
input symmetry. The mandatory full-alphabet dual comparison and the
suspension verdict above are unchanged.

## Completed correction — 2026-08-31: exact orbit-total replay restores all 20 $q=3$ rows

Campaign:
`campaigns/2026-08-31T09:13:25Z_7dc5babe-5e1b-44fa-a9e5-03d2f10183b6_q3-invariance-correction/`;
manifest: `manifest.json`. Frozen correction runner SHA-256
`ed84c1a3982f243e49e8f8f35367ee9e8788bb3c8a47fed16b883bd9da9efa21`;
frozen conservative renderer SHA-256
`bc8901f1038d09650fae9dea1caaaab25861f5102ae595efa3e17601782e1b15`.

**MACHINE-VERIFIED correction verdict:** exact audit found that 16/20 frozen
accepted `ba_word` vectors split an output orbit; those 16
representative-compressed upper endpoints remain retracted. The other four
legacy vectors were invariant. The replacement snaps total input-orbit masses
for the primal and total output-orbit masses for the dual, distributes each
mass uniformly inside its orbit, proves exact invariance, and compares both
the primal and dual orbit evaluations against direct evaluation over every
word of the full alphabet. Every chosen dual was
`ba_total_orbit_mass`. All 20 replacements certify
`CERT_LOWER_BEATS_LBplus + CERT_UPPER_BEATS_UB`; no
Tavakoli–Nguyen–Bose theorem or printed sandwich value is refuted.

The table below uses only conservative outward decimals derived from the full
400-bit per-symbol Arb balls. The binary64 `cert_lo_per_symbol`,
`cert_hi_per_symbol`, and width fields in the detailed JSONL are
**REPORT-ONLY**, not certificate endpoints.

| $n$ | $d$ | corrected certified interval (outward) | width upper (outward) | legacy $D$ split? |
|---:|:---:|:---|---:|:---:|
| 6 | 1/2 | `[0.48612690873525893, 0.48612726836569464]` | `3.5963043559742424e-07` | no |
| 6 | 1/5 | `[0.9992805994038761, 0.9992806185950381]` | `1.919116167554878e-08` | yes |
| 6 | 1/10 | `[1.262451538526818, 1.2624515509732321]` | `1.2446413657816448e-08` | yes |
| 6 | 1/20 | `[1.4158072326903546, 1.4158072449168306]` | `1.2226475645090283e-08` | yes |
| 7 | 1/2 | `[0.46572596820424017, 0.46573204857053757]` | `6.080366297333852e-06` | yes |
| 7 | 1/5 | `[0.9769505177740003, 0.976950563790745]` | `4.6016744441483404e-08` | yes |
| 7 | 1/10 | `[1.2473597256274749, 1.247359775851469]` | `5.022399355970927e-08` | yes |
| 7 | 1/20 | `[1.4070914100656966, 1.4070914501217517]` | `4.005605471757835e-08` | no |
| 8 | 1/2 | `[0.4491602653164887, 0.44918246544877494]` | `2.220013228610702e-05` | yes |
| 8 | 1/5 | `[0.9582453246473193, 0.9582455337027981]` | `2.09055478620426e-07` | yes |
| 8 | 1/10 | `[1.23445208950039, 1.2344521805124733]` | `9.101208270304075e-08` | yes |
| 8 | 1/20 | `[1.399541506237329, 1.399541571855359]` | `6.56180294580417e-08` | yes |
| 9 | 1/2 | `[0.43538141101923905, 0.4353817875980781]` | `3.765788389715567e-07` | yes |
| 9 | 1/5 | `[0.9423194816428748, 0.9423200290509085]` | `5.474080334643068e-07` | no |
| 9 | 1/10 | `[1.2232499579738394, 1.2232502176397095]` | `2.5966586972603587e-07` | yes |
| 9 | 1/20 | `[1.3929068203510147, 1.3929070242522732]` | `2.03901257961709e-07` | yes |
| 10 | 1/2 | `[0.4236997564364154, 0.42370155403328696]` | `1.7975968714839528e-06` | no |
| 10 | 1/5 | `[0.9285797540357409, 0.9285822403778092]` | `2.4863420681301732e-06` | yes |
| 10 | 1/10 | `[1.2134129645213867, 1.2134136998846528]` | `7.353632655671914e-07` | yes |
| 10 | 1/20 | `[1.3870088535516931, 1.387009425126998]` | `5.715753042954663e-07` | yes |

The replay performed 1,058,508 direct full-input KL evaluations, 5,821,704
exact group-generator checks, and 237,744,816 positive-$W$/positive-$D$ term
checks (`MACHINE-VERIFIED` counts); zero positive-mass input words were
skipped. Runtime `841.039 s` wall / `838.785 s` CPU is
`COMPUTATIONAL-EVIDENCE`.

**Rule-7 scope sentence.** The correction sweep is exactly $q=3$ with
$(n,d)$ in
`{6,7,8,9,10} x {1/2,1/5,1/10,1/20}`, using
$G=S_3\times C_2$ on symbol VALUES and reversal, exact total input-orbit mass
snap $2^{30}$, exact total output-orbit mass snap $2^{30}$, full finite
channel support for every dual candidate, direct full input alphabets of size
$3^n$, and Arb 400-bit outward rounding. **NOT swept:** $q=4$ or $q\ge4$;
$q=3,n\ge11$; other $d$; per-word tie-breaking as an admitted bound;
non-$G$-invariant correction candidates; the cause of Tavakoli–Nguyen–Bose
rounding; Morozov–Duman, Pinto–Ribeiro, and every $n\to\infty$/asymptotic
capacity claim. Every replacement is only a finite-$n$ theorem.

Disproved premises: equal orbit-invariant float inputs do **not**
automatically remain equal under stable per-word largest-remainder snapping,
and only 4—not 20—of the frozen accepted `ba_word` vectors were invariant.
The hazard was found retrospectively in the parent $q=4$ startup audit; the
separate correction box, orbit-total formulas, exact checks, candidate order,
and endpoint rule were prospectively frozen before its runner.

The parent $q=4$ campaign remains **FAILURE TO CERTIFY / NOT RUN**: zero
$q=4$ certificate rows executed, pending owner integration. Named next
campaign: **q4-total-output-orbit-mass restart**, replaying $q=4,n=5..10$ at
the four fixed $d$ values and then separately adjudicating $q=3,n=11$.
**INFERENCE cost:** 6–10 single-thread CPU hours under the existing 96 GiB,
90-minute-per-row, and 10-hour-stage caps.

## $q=4$ total-output-orbit-mass restart — static freeze, 2026-08-31

Campaign:
`campaigns/2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/`.
The exact 24-row box is preregistered and runnable, but remains **NOT RUN**:
zero BA iterations, zero deferred-anchor rows, and zero $q=4$ certificate rows
were executed while another target owned the heavy CPU slot
(`COMPUTATIONAL-EVIDENCE`). The frozen authoritative runner SHA-256 is
`469e31a7cb59712d0eea6e5d932681fce63ff8eecd78b093d2da37be89ae09d7`;
the static checksum-ledger SHA-256 is
`906702e1d649cda653a7058792d274fe7b3cbdc88b830f793daa5a6d80aa999d`
(`MACHINE-VERIFIED` hashes).

The construction admits only the corrected input-orbit-total primal and
output-orbit-total dual candidates `(ba_total_orbit_mass, uniform)`;
`ba_word` is forbidden. Every future row requires exact sum,
positivity/support and $S_4\times C_2$ invariance checks, direct
full-alphabet KL comparison over all $4^n$ inputs, 400-bit Arb evaluation,
and exact containment by conservative outward decimals.

The sub-second static guard passed (`MACHINE-VERIFIED`): the $q=4,n=5$
census tuple is `(31,50,631,332,963)`, 130,560 value/reversal equivariance
checks passed, the position plant broke with $A=1$ versus $A=2$, the planted
split output law was rejected before KL, zero support returned infinity, and
the negative-width plant was rejected. Empty resume validation found exactly
zero rows.

Provenance correction: the row JSONL and checksum-ledger replacements are
each independently atomic and `fsync`ed, but the **pair is not
crash-atomic**. A crash between them preserves a mismatch and hard-aborts
resume; it is never accepted or repaired. Failed temporary evidence is never
unlinked and its exact path is reported.

**Rule-7 scope sentence.** The intended sweep is exactly $q=4$ with
$(n,d)$ in
`{5,6,7,8,9,10} x {1/2,1/5,1/10,1/20}`, using
$G=S_4\times C_2$ on symbol VALUES and reversal, total input-orbit and
output-orbit snaps $2^{30}$, the two fixed dual candidates above, direct
full-alphabet KL, full support, and Arb 400-bit outward rounding. **NOT
swept:** $q=3,n=11$ or any q=3 row; $q\ge5$; $q=4,n\ge11$; other $d$;
`ba_word`; non-invariant candidates; the falsified position/type reduction;
the cause of Tavakoli–Nguyen–Bose rounding; Morozov–Duman, Pinto–Ribeiro, and
every asymptotic-capacity claim. The deferred next action is the frozen
$q=4,n=5,d=1/2$ full startup anchor after Main releases the heavy slot,
followed only on PASS by the fixed 24-row run. **INFERENCE cost:** 6–10
single-thread CPU hours.

## $q=4$ released anchor passes; 24-row grid still pending (owner, 2026-08-31)

The registered full startup anchor in
`campaigns/2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/`
has now run under the final source hash `c9c3a842…`. A fresh release-time
static guard passed first. The exact census again matched
`(31,50,631,332,963)`; input/output orbit sizes sum to 1,024/1,365.

For `q=4,n=5,d=1/2`, the admitted exact total-orbit primal and
`ba_total_orbit_mass` dual give the outward per-symbol interval

`[0.6664806108007938, 0.6664806978587714]`,

width upper `8.705797740570955e-8 < 1/500`. All 3,072 input generator checks
and 4,095 checks for each output candidate passed, and both dual candidates
were evaluated on all 1,024 inputs. Verdict:
`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`
(**MACHINE-VERIFIED**).

Anchor SHA-256:
`c47aeb35250eade511e3cfe8e65701e9302545b18d46d995952cdd5f073b794b`.
The anchor deliberately wrote no row artifact. Therefore this restores only
the registered anchor: the full 24-row
`n=5..10`, `d in {1/2,1/5,1/10,1/20}` production grid remains
**COMPUTE PENDING**, and no asymptotic claim changes.

## $q=4$ 24-row grid completed — exact orbit-total certificates (owner, 2026-08-31)

This append-only result supersedes the immediately preceding
**COMPUTE PENDING** status. The frozen `c9c3a842…` runner completed all 24
registered rows in
`campaigns/2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/`.
Every row is `CERTIFIED`; every selected dual is
`ba_total_orbit_mass`; and every conservative interval strictly beats both
published finite-$n$ endpoints
(`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`,
**MACHINE-VERIFIED**).

The complete interval table is
`TABLE_q4_conservative_intervals.csv`. Width upper bounds range from
`6.165071504436472e-09` at $(n,d)=(5,1/10)$ through
`2.3125135592425783e-05` at $(8,1/2)$ bits/symbol. The tightest published
margin remains positive at $(5,1/20)$: the serialized Arb lower-margin ball
has lower endpoint approximately `0.0007402311927961205436222872684`, and
the upper-margin ball has lower endpoint approximately
`0.08766315500354434284029482356`. The archived exact Arb balls, not these
decimal summaries, certify the strict inequalities.

Across the fixed pair of dual candidates, the rows record 11,182,080 direct
full-input evaluations and 3,585,643,216 positive-$W$/positive-$D$ term
checks. The primal side records 1,792,821,608 direct conditional entries.
All 16,773,120 expanded input and 44,728,272 expanded output generator
checks passed; there were zero full-word/representative overlap failures.
Post-run immutable-resume validation returned `PASS` with 24 rows. An
independent owner check verified all 24 newline-inclusive row hashes and
that the CSV is an exact projection of the canonical JSONL
(**MACHINE-VERIFIED**).

The one-process production run took `4651.466 s` wall / `4649.094 s` CPU;
maximum row-recorded RSS was `1,441,169,408` bytes and no resource stop
fired (**COMPUTATIONAL-EVIDENCE**). Full report:
`report.md`; additive final manifest:
`manifest_final.json`; final checksum-ledger SHA-256:
`d2660432da9a4d79fb3d9c0e134a249cef6c546fe199f2a3bce672dce905eb30`.
The historical manifest, interim ledger, and post-hold provenance amendments
remain byte-preserved; import-generated `__pycache__/` files are preserved
but excluded from evidence ledgers.

**Rule-7 scope sentence.** This proves only the registered finite-$n$ box
$q=4$, $n=5..10$, and
$d\in\{1/2,1/5,1/10,1/20\}$ under the frozen
$S_4\times C_2$ exact-orbit protocol. It does **not** sweep $q=3,n=11$,
$q\ge5$, $q=4,n\ge11$, other deletion probabilities, non-invariant
candidates, or any asymptotic-capacity statement. Named next campaign:
**$q=3,n=11$ total-output-orbit-mass extension**, with a separate
prospective freeze required before computation.

## $q=3,n=11$ extension completed — 4/4 orbit-total certificates (owner, 2026-08-31)

The named successor campaign is complete. The frozen
`08ecbc49…` runner certified all four registered rows in
`campaigns/2026-08-31T20-34-52Z_04ed25bd-455f-4814-8288-5ee0106e0db8_q3-n11-total-output-orbit-mass/`.
Every row is `CERTIFIED`, every selected dual is
`ba_total_orbit_mass`, and every conservative interval strictly beats both
published Tavakoli–Nguyen–Bose finite-$n$ endpoints
(`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`, **MACHINE-VERIFIED**):

| $d$ | certified interval (outward) | width upper |
|:---:|:---|---:|
| 1/2 | `[0.4136415126778377, 0.4136475665528263]` | `6.053874988461827e-06` |
| 1/5 | `[0.9165946848818924, 0.9166017532024916]` | `7.068320598920295e-06` |
| 1/10 | `[1.2046904010391062, 1.2046925848595933]` | `2.1838204867062363e-06` |
| 1/20 | `[1.3817159992481873, 1.3817177042129642]` | `1.7049647766059455e-06` |

All four rows reproduce the exact census tuple
`(14884, 22450, 6148309, 3573542, 9721851)` with orbit-size sums
$177147=3^{11}$ and $265720=\sum_{k=0}^{11}3^k$, matching independent
Burnside counts. The rows record 1,417,176 direct full-input dual
evaluations, 634,415,384 positive-$W$/positive-$D$ dual terms,
316,716,419 direct primal conditional entries, 2,125,764 expanded input
generator checks, and 9,565,920 expanded output generator checks, with zero
full-word/representative overlap failures. The frozen full-support bump
fired on `ba_total_orbit_mass` at $d=1/2,1/10,1/20$
(5,448/1/4 zero orbit masses) and not at $d=1/5$. The tightest published
margin is at $d=1/20$, where the archived Arb lower-margin ball has lower
endpoint approximately `0.00129169637501813941901314966788`.

Post-run immutable resume returned `PASS` with four rows; an independent
owner check re-derived all four newline-inclusive row hashes, re-verified
every outward decimal against its archived exact binary rational, re-checked
`width <= 1/500`, and confirmed the CSV is an exact projection
(**MACHINE-VERIFIED**). The one-process run took `1151.072 s` wall /
`1150.401 s` CPU with maximum row-recorded RSS `929,677,312` bytes and no
resource stop (**COMPUTATIONAL-EVIDENCE**). Report: `report.md`; additive
`manifest_final.json`; final checksum-ledger SHA-256
`0ba5629a13854ba4eb1cdff32aad83fa0a385ec9a3da6b33bef5f94ac28b7761`, all 20
entries passing. No `__pycache__` exists: bytecode writing is refused before
the first campaign import.

**Rule-7 scope sentence.** This proves only the registered finite-$n$ box
$q=3$, $n=11$,
$d\in\{1/2,1/5,1/10,1/20\}$ under the frozen
$S_3\times C_2$ exact-orbit protocol. It does **not** sweep any $q\ge4$ row,
$q=3$ with $n\ne11$, other deletion probabilities, `ba_word`, non-invariant
candidates, or any asymptotic-capacity statement. Named next campaign:
**$q=3,n=12$ total-output-orbit-mass extension** (input orbits 44,530 over
$3^{12}=531441$ words; output orbits 66,980 over
$\sum_{k=0}^{12}3^k=797161$ words), requiring its own prospective freeze.
**INFERENCE cost:** about 1.0–1.5 single-thread CPU hours for its four rows
under the same caps.

## Current state (agent DelcapN12, 2026-09-02) — $q=3,n=12$ box COMPLETE, 4/4 certified

The $q=3,n=12$ extension named above is **done**: all four registered rows
are `CERTIFIED` and **every one strictly beats BOTH published
Tavakoli–Nguyen–Bose finite-$n$ endpoints**
(`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`, `MACHINE-VERIFIED`). Finite-$n$
theorem for exactly those four cells; **no asymptotic-capacity claim.**

Two campaigns, because the first crashed on a launcher defect (below):

* `campaigns/20260901T121346Z_d1db3155_4354a957b8ba/` — gate
  `delcap-q3-n12-total-output-orbit-mass-v1`, verdict **CRASHED**. Box rows 0–1
  ($d=1/2,1/5$) certified and ledger-validated before the crash.
* `campaigns/20260902T010546Z_adc6518e_e85083dd9888/` — gate
  `delcap-q3-n12-continuation-d10-d20-v1`, verdict **FROZEN-CERTIFIED**. Box
  rows 2–3 ($d=1/10,1/20$). Canonical report: that dir's `report.md`.

| box row | $d$ | certified interval (outward) | width upper | $LB^+$ margin $\ge$ | $UB$ margin $\ge$ |
|---:|:---:|:---|---:|---:|---:|
| 0 | 1/2 | `[0.40486950906723745, 0.4048828292730583]` | `1.3320205820702559e-05` | `9.774185e-02` | `3.875984e-01` |
| 1 | 1/5 | `[0.9060415197452321, 0.9060644026068896]` | `2.2882861657318862e-05` | `2.095869e-02` | `3.619056e-01` |
| 2 | 1/10 | `[1.1968923784214915, 1.1968988416160002]` | `6.4631945083381036e-06` | `5.290779e-03` | `2.295674e-01` |
| 3 | 1/20 | `[1.376928469567818, 1.3769335302038865]` | `5.060636068201896e-06` | `1.314571e-03` | `1.287808e-01` |

Every selected dual is the exact `ba_total_orbit_mass` candidate; exact
orbit-uniform was evaluated over the full input alphabet in every row and never
selected. All four widths are inside the registered $1/500$ target.

**Pre-registered orbit counts asserted, never adjusted** (`MACHINE-VERIFIED`):
input $S_3\times C_2$ orbits $44{,}530$ over $3^{12}=531{,}441$ words; output
orbits $66{,}980$ over $\sum_{k=0}^{12}3^k=797{,}161$; independent Burnside
closed form $B(k)=(3^k+3+\mathrm{rev}(k))/12$ reproduces both, with series
$B(0..12)=1,1,2,4,10,25,70,196,574,1681,5002,14884,44530$ whose $k\le11$ prefix
is $22{,}450$ — the frozen $n=11$ anchor. Full census tuple
$(44530, 66980, 30666848, 18311678, 48978526)$, identical across both static
guards and all four rows. **Exact orbit-mass identities:** input numerators sum
to exactly $2^{30}$; the BA dual sums to $2^{30}$ with no bump ($d=1/5$) and to
$2^{30}+66{,}980=1{,}073{,}808{,}804$ where the frozen full-support bump fires
($30904$, $2$, $8$ zero masses at $d=1/2,1/10,1/20$); the orbit-uniform dual's
numerators sum to the total output word count $797{,}161$.

**Controls, in-run, both directions.** ACCEPT: the frozen $n=11$ row
`0:3:11:1/2` was re-executed through the $n=12$ code path and reproduced
**bit-exactly on all 35 compared certificate fields** — census tuple, all
$14{,}884$ input and $22{,}450$ output mass numerators, every archived Arb ball,
the whole `conservative_interval` block with its exact binary rationals, the
verdict, and even the float BA locator values — zero mismatches, in both
campaigns. REJECT: five plants each fired at their pre-registered assertion
point (`R1a` mass total $2^{30}+1$; `R1c` negative input numerator; `R1b`
negative output numerator; `R2` dual-below-primal negative width; `R3` planted
split non-invariant output law rejected before compression with 5 failing
generator equalities, while the exact orbit-uniform control was admitted).

Aggregate exact work over the four rows: $1{,}553{,}383{,}596$ direct primal
conditional entries, $4{,}251{,}528$ direct full-alphabet dual word evaluations,
$6{,}377{,}292$ expanded input generator checks, $28{,}697{,}796$ expanded output
generator checks, zero overlap failures, $32{,}442$ zero-mass input words
accounted for at $d=1/2$. Row CPU $7664.5$ s $=2.13$ CPU-hours inside the
pre-registered 6-CPU-hour budget; max row RSS $4{,}398{,}972{,}928$ bytes against
the 96 GiB cap; no resource stop (`COMPUTATIONAL-EVIDENCE`). A 180-assertion
independent re-derivation from the frozen bytes of both campaigns passed at
close-out (line hashes, canonical-JSON equality, outward-rational ↔ decimal
identity, outward containment against the archived Arb balls, width $\le 1/500$,
both-endpoint strict improvement, census determinism, mass identities).

**Instrument defect disclosed.** The first campaign died mid-row-2 with
`BrokenPipeError(32)` raised inside the runner's `log()` `print`: the launcher
piped stdout through `tail`, deviating from the frozen command, and the read end
closed. Root cause is the launcher, not the mathematics; the runner's frozen
semantics correctly hard-aborted, archived `FAILURE.json`, and forbade rerun, so
that campaign was frozen and closed `CRASHED` with its two valid rows intact and
byte-pinned. The residual defect worth carrying forward: `log()` is fatal on a
broken stdout even though all evidence is already durable on disk. The
continuation excludes the mode by contract (stdout to a file, supervised
process) rather than by patching frozen mathematics. Also disclosed: two
close-out *checking* errors of mine (comparing the archived outward width to the
difference of display doubles, ~2 ULP larger than the true Arb difference the
width bounds; and expecting the post-bump denominator to be $2^{30}+z_{\text{zeros}}$
rather than $2^{30}+n_{\text{out}}$) — both were wrong checks, not wrong evidence.

**Next gate:** $q=3,n=13$ total-output-orbit-mass extension under its own
prospective freeze — input orbits $B(13)=133{,}225$ over $3^{13}=1{,}594{,}323$
words, output orbits $\sum_{k=0}^{13}B(k)=200{,}205$ over
$\sum_{k=0}^{13}3^k=2{,}391{,}484$ words. **INFERENCE cost:** the measured
$n{=}11\to n{=}12$ row-CPU step was $1150\,\mathrm{s}\to 7664\,\mathrm{s}$
(6.7x), so $n=13$ needs about $6$–$12$ single-thread CPU-hours for four rows and
roughly $12$–$14$ GiB peak RSS; it must be pre-registered with a larger budget,
never resumed under changed caps.

## Current state (agent DelcapN12, 2026-09-02b) — $q=3,n=13$ box COMPLETE, 4/4 certified

The $q=3,n=13$ extension named above is **done**: all four registered rows are
`CERTIFIED` and **every one strictly beats BOTH published Tavakoli–Nguyen–Bose
finite-$n$ endpoints** (`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`,
`MACHINE-VERIFIED`). Finite-$n$ theorem for exactly those four cells; **no
asymptotic-capacity claim.** One campaign, single clean run, no crash, budget
not binding:

* `campaigns/20260902T024656Z_90c757de_c4e9905d05a1/` — gate
  `delcap-q3-n13-total-output-orbit-mass-v1`, verdict **FROZEN-CERTIFIED**.
  Canonical report: that dir's `report.md`.

| $d$ | certified interval (outward) | width upper | $LB^+$ margin $\ge$ | $UB$ margin $\ge$ |
|:---:|:---|---:|---:|---:|
| 1/2 | `[0.39713666568100003, 0.397161636764393]` | `2.4971083392840876e-05` | `9.859301e-02` | `3.953196e-01` |
| 1/5 | `[0.8966732662961678, 0.8967669676702715]` | `9.370137410351712e-05` | `2.129289e-02` | `3.712030e-01` |
| 1/10 | `[1.1898717660767768, 1.1898901620041096]` | `1.8395927332359186e-05` | `5.372100e-03` | `2.365761e-01` |
| 1/20 | `[1.3725688358567392, 1.3725840117520043]` | `1.5175895264739525e-05` | `1.333958e-03` | `1.331304e-01` |

Every selected dual is the exact `ba_total_orbit_mass` candidate; exact
orbit-uniform was evaluated over the full input alphabet in every row and never
selected. All four widths inside the registered $1/500$ target.

**Pre-registered orbit counts asserted, never adjusted** (`MACHINE-VERIFIED`):
input orbits $B(13)=133{,}225$ over $3^{13}=1{,}594{,}323$ words; output orbits
$\sum_{k=0}^{13}B(k)=200{,}205$ over $\sum_{k=0}^{13}3^k=2{,}391{,}484$;
Burnside series $B(0..13)$ reproduces both, and its $k\le11$/$k\le12$ prefixes
($22{,}450$/$66{,}980$) reproduce the frozen $n=11$/$n=12$ anchors. Measured
dense census tuple $(133225, 200205, 152963378, 93991184, 246954562)$,
byte-identical across the static guard, both production structure builds, and
all four rows; the pre-registered projection ($\sim$1.53e8/9.4e7/2.5e8) matched
to $\le$0.2%. **Exact orbit-mass identities:** input numerators sum to $2^{30}$;
BA dual sums to $2^{30}$ (no bump, $d=1/5$) or $2^{30}+200{,}205=1{,}073{,}942{,}029$
where the full-support bump fires ($131{,}380$/$4$/$18$ zeros at
$d=1/2,1/10,1/20$); orbit-uniform numerators sum to the total output word count
$2{,}391{,}484$.

**Controls, in-run, both directions:** ACCEPT — the frozen $n=12$ row
`0:3:12:1/10` re-executed through the $n=13$ code path reproduced **bit-exactly
on all 41 compared certificate fields** (census tuple, all $44{,}530$ input and
$66{,}980$ output numerators, every Arb ball, the whole conservative-interval
block with exact binary rationals, the verdict, even the float BA locator
$14.362708541093044$), zero mismatches; its line hash was ledger-verified
first. REJECT — R1a ($2^{30}+1$ input total), R1b (negative output numerator),
R1c (negative input numerator), R2 (negative width), R3 (planted split output
law, rejected before compression with failing generator equalities while the
orbit-uniform control was admitted) — each fired at its pre-registered
assertion point.

Aggregate exact work: $7{,}446{,}791{,}969$ direct primal conditional entries,
$19{,}131{,}876$ expanded input and $86{,}093{,}424$ expanded output generator
checks, zero overlap failures, $289{,}980$ zero-mass input words accounted at
$d=1/2$. Stage $38{,}668.4$ s CPU = **10.74 of the registered 24 CPU-hours**,
wall 11.28 h of 26, max row RSS $16{,}669{,}851{,}648$ bytes of the 32 GiB cap
(48.5%); static-guard build peak $16.83$ GB (49.0%) — **the budget did not
bind and no resource stop occurred.** The pre-registered launcher rule (no
stdout piping; the $n=12$ BrokenPipe lesson) was obeyed and the failure mode
never recurred: single clean run, exit 0, no `FAILURE.json`. A 180-assertion
independent re-derivation from the frozen bytes passed at close-out; both
freeze ledgers re-derive with zero mismatches. A per-function AST diff proves
22 certificate/ledger functions byte-identical to the $n=12$ instrument
(changed only: `rows_stage` row-count assertion and the census-target
constants; new: the $n=12$ ACCEPT loader/comparator and the census-pin
helper).

**Next gate: $q=3,n=14$ — evaluated, NOT opened.** Measured row-CPU scaling
$287.6$ s ($n{=}11$) $\to$ $1{,}916.1$ s ($n{=}12$) $\to$ $9{,}573.8$ s
($n{=}13$) projects $n=14$ at $\sim$$13.9$ h/row, $\sim$$55.5$ CPU-hours for
four rows, multi-day stage wall, and a $\sim$$35$–$40$ GB live-structure peak
(measured $n{=}13$ RSS $16.67$ GB at $2.47\times10^8$ sparse slots) — above
the 32 GiB envelope this workstation campaign family pre-registers, with real
thrash risk against sibling campaigns. Structural targets if it is ever
opened: $B(14)=399{,}310$ input orbits over $3^{14}=4{,}782{,}969$ words,
$\sum_{k=0}^{14}B(k)=599{,}515$ output orbits over
$\sum_{k=0}^{14}3^k=7{,}174{,}453$ (machine-recomputed; an earlier draft of
this line said $394{,}024$/$594{,}229$, wrong, corrected). Opening $n=14$
requires a fresh preregistration with a $\sim$64 GiB RSS envelope (probe and
resource-stop rule retained), a $\sim$64 CPU-h budget, and a stated position
on the multi-day single-thread wall — or a separately frozen instrument
change (chunked/array-backed structure build). Per the pre-registered rule,
the projection is recorded and the campaign family stops here at a complete
certified box.
