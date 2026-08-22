# Hilbert's Tenth Problem over Q (`h10q/`)

**Target.** Is the existential (positive) theory of $\mathbb{Q}$ decidable —
equivalently, is there an algorithm deciding, given $f_1,\dots,f_m \in
\mathbb{Q}[x_1,\dots,x_n]$, whether they have a common rational zero? Open in
both directions. The 2024–25 breakthroughs (Koymans–Pagano arXiv:2412.01768;
Alpöge–Bhargava–Ho–Shnidman arXiv:2501.18774) settled every ring of integers;
the frontier is $\mathbb{Q}$ itself (Koymans–Pagano expository arXiv:2602.04468,
Question 6.2). Opened 2026-08-12 in `scratch/`, canonicalized here 2026-08-15;
nothing in `scratch/` remains.

## Authorities (single source of truth)

| File | Owns |
|---|---|
| [`THEOREMS.md`](THEOREMS.md) | Precise frontier statements (L1–L9, A1–A3), full proofs of everything proved here, and the exact open assembly lemma. |
| [`RESULTS.md`](RESULTS.md) | Findings ledger, claim-by-claim scope table, frontier stratification, next actions. |
| [`NOTES.md`](NOTES.md) | Problem map: exact theorem statements with arXiv ids, both pillars. |
| [`h10q.py`](h10q.py) | Standalone executable model; re-verifies every in-suite claim on each run. |
| [`l6_search.py`](l6_search.py) | Guarded tied-certificate searcher (`--canonical` enforces the W0/W1 hypotheses). |
| [`l9_steer.py`](l9_steer.py) | Reciprocity-steered assembly searcher + deep-rescue mode; certificate logic lives in `h10q._l9_steered_solvable` (single source). |
| [`l13_filter.py`](l13_filter.py) | Zero-bad square/smooth filter (`emergent_free`: strip the frozen set ∪ supp(b) from P(b), test remainders are squares — no factorization, no engine refusals) + the `cofactor_decide` ladder (square / proved-prime / prime-power / factorint rungs + Jacobi evidence tier) + cell-(89,(-2,1)) attack driver. Search-side: exact counters only, never citable as cell coverage. |
| [`l13h_scan.py`](l13h_scan.py) | Per-class emergent-free member certification for hypothesis H (all 293 grid classes; `smooth_emergent` + `cofactor_decide`, stop-on-first-zero, reason-coded refusals). Data: `data/l13h_{l11_a,l11_b,esc}[_2].jsonl`. |
| [`test_bls.py`](test_bls.py) | Manual unit check of the BLS-relaxed Pocklington (pinned-regime proofs + composite-cousin guard; NOT part of the suite). |
| [`l13h_alt.py`](l13h_alt.py) | Wave-3 alternate-class driver: Hensel-root a-lifts (L11), f/ε variants (ESC, self-validating constructor), DEEPK mode; a=1 fallback for w≡3 mod 4. Lead-owned. |
| `data/l13h_{l11_a,l11_b,esc}{,_2}.jsonl` | H-scan waves 1–2 per-class records (frozen classes, k≤200): verdict histograms + closures, reason-coded refusals. |
| `data/l13h_alt_{l11_a,l11_b,esc,l11_hard,esc_hard}.jsonl`, `data/l13h_deepk_l11.jsonl` | Wave-3 alternate-class + deep-k closure records. |
| `data/l13h_qa.{jsonl,md}` | LadderQA regression: 45/45 exact agreement (29 zero + 16 bad rows), 0 alarms. |
| `data/l13h_all_closures.json` + `data/l13h_replay.jsonl` | Canonical merged 293-cell closure set + lead-wide replay verdicts (293/293 ladder-zero; ramified 271 empty / 22 refused-logged). Authority behind `_L13_H_CLASSES`. |
| `data/l15_remainders.jsonl` | L15 exact remainder audit over all 293 closure members (R squarefree in every row; 197 $R=p$ exactly, 96 factorint $e\in\{2,3,4\}$; 0 refusals, 0 alarms). |
| `data/l15_density.json` | L15 first-member counting law (deterministic 2.8-s playback): $k_{\text{zero}}$ histogram, per-family/per-rung splits, Bateman–Horn calibration, outlier table. |
| `data/l16_char.jsonl` | L16 probe (lead-run): full sign-decorated emergent lists per prime-member $k$ on 7 sampled classes; residue-determination (small layer sieve-shaped) + two-layer obstruction split; parity row-checks. |
| `data/l17_sieve.jsonl` | L17: all-293-class member-level obstruction scan ($k\le119$, 8,760 prime members): residue-determination 293/293, clean-rate histograms, recurring-prime tables. Persisted evidence (see **Reproducibility scope** below). |
| `data/l17_badroots.jsonl`, `data/l17_badroots_closures.jsonl` | L17: exact bad-signed root residues ($k\bmod p$) with certified odd-valuation masses — wave-1 cohort and closure-authority versions; closure-authority reconciled 5,375/5,375 against observed occurrences. (The wave-1 cohort variant `l17_badroots.jsonl` is persisted evidence without a checked-in producer.) |
| [`l17_ratemodel_padic.py`](l17_ratemodel_padic.py) + [`l17_badroots_closures.py`](l17_badroots_closures.py) + [`l17_matched.py`](l17_matched.py) + `data/l17_{badroots_closures,matched}.jsonl` | L17: matched-test regenerator + statistics: marginal ratio 0.9993; union residual −2.9/8760; per-class rows; marginal-preserving permutation $z=-0.09$, empirical tails $0.61/0.53$. Class-factors and the Haar-window estimand are $p\le10^4$-window statistics only (no $p>10^4$ tail bound; the excluded-mass residual covers singular root lifts). |
| `data/l17_ratemodel_censored.jsonl` | L17: two-layer rate model on the predeclared wave-1 cohort, right-censored with two-sided $c$-bands (unknown trials counted both ways); family $c \in [0.017, 0.80]$. Persisted evidence; no checked-in producer. |
| `data/l17_cofactor.jsonl` | L17: cofactor sign audit on closures: proved reciprocity-level statements (197/197 auto-$+1$; parity product $+1$ on 293/293); factorint/residue tables closure-conditioned. Persisted evidence; no checked-in producer. |
| `data/l17_cofactor_nonclosure.jsonl` | L17: nonclosure cofactor sign sample, power-gated: 11/400 resolved (389 prover refusals); parity-consistent; distributional claim INCONCLUSIVE (declared under-powered in the global summary). Persisted evidence; no checked-in producer. |
| `data/l17_classfactors.jsonl` | L17: per-class $p\le10^4$ window factor intervals [lo, hi] rebuilt from v2 resolved-mass roots (generator [`l17_classrisk.py`](l17_classrisk.py)); closure-effort correlation (descriptive): Spearman $\rho=-0.178$, medians $0.5742$ vs $0.5662$ (cls=0 vs ≥1). |
| `data/l17_schinzel_audit.jsonl` | L17: audit of the ledgered Schinzel reduction (THEOREMS L9–L10): 293 canonical closure pairs $\{q_1(t), F_{\text{cell}}(t)\}$ — all $F$ irreducible (deg 8, certified), value-gcd 1, pair-product gcd 1, local conditions pass with corrected $q_1(t)\not\equiv0\ (p)$ counts; 6 `_L13_RESIDUAL` rational-b witnesses segregated as pointwise certificates. Generator [`l17_schinzel_audit.py`](l17_schinzel_audit.py). |
| `data/l17_horizon{101,103,107,109}.jsonl` + [`l17_horizon{101,103,107,109}.py`](l17_horizon101.py) | L17-horizon: off-grid closure artifacts + full scan logs at $w\in\{101,103,107,109\}$, $u=-1$; replay from repo root; H109 cross-k auxiliary scan `data/l17_horizon109_crossk.jsonl` (1,275 probes, $k\le50$; deterministic generator [`l17_horizon109_crossk.py`](l17_horizon109_crossk.py), ~17 min, byte-identical regeneration verified). |
| Engine note (2026-08-18) | Proven-primality engine strengthened: `_pocklington` now uses trial division to $10^6$ (early abort at $F\ge n^{1/3}$, exact integer cube root) + the **Brillhart–Lehmer–Selfridge relaxation** (CP Thm 4.1.5: $F\ge n^{1/3}$, two-factor discriminant test). Every verdict still exact; `PrimalityBound` refusals remain non-evidence. Suites: default ~51 s, extended ~118 s. |
| `data/l6_witnesses.jsonl`, `data/l9_steered.jsonl` | **Canonical serializations** of the authorities `_L6_WITNESSES` (61 rows) and `_L9_STEERED` (345 rows) in `h10q.py`. Byte-identical match required on every suite run — drift, absence, or any underived field fails the suite. Regenerate: `python3 h10q.py --export-evidence`. |
| `data/l9_steer_run.jsonl`, `data/l9_rescue_*.jsonl` | Raw steered-search run artifacts — search-side provenance only, never citable evidence. |
| `data/l13_filter_run2.json` | L13f stage-2 outcome (2026-08-18): full box spec, 22 counters with reason-coded refusals, 4,481 stored cofactor rows with per-row verdicts, 236 hit records, dense_box block. Authority for the cell-closure counters; the suite-frozen witness is the closing claim. |
| `l13_dense.py` + `data/l13_dense_smallb.jsonl` + `data/l13_dense_cofactorbig.jsonl` | DenseScan dense deep box: every odd $|b|<20001$, four branches, $a\in\{1,17,25,33\}$ — 319,968 candidates, 3,041 aligned (ALL $a{=}17$), 1,286 cofactor-big rows (handoff to the L13f ladder → 89 certified soluble). Exact counters only. |
| `data/litscout_h10q.md` | LitScout-2 (2026-08-18): source-verified citations for CONDITIONAL.md §4 — Krumm JTNB 28 (2016) pins (Thm 1.3, Props 3.4–3.5; degree 3 conditional on the elliptic Parity Conjecture, Prop 3.8), degree-8 ceiling beyond proved range on both sides, DDF arXiv:2102.06941v5 attribution confirmed, Crelle 495 (1998) author-list correction. Applied to CONDITIONAL.md same day. |
| `data/superseded/` | Inadmissible or superseded raw search buckets and the old rescue log — provenance only, never citable evidence. |

## State (2026-08-17)

- **PROVED here (machine-checked):** W0 canonical two-branch selector from the
  character identity $\sum_{a\in k}\chi(1+4a^2)=-1$ over every odd residue
  field including $\mathbb{F}_3$; W1 exact tied target-place criterion
  $\chi_w(-\delta_\tau A)=1$; W2 exact global criterion — the disjunction
  $M_\tau=0$ **or** $-\delta_\tau A M_\tau\in N_{E_\tau/\mathbb{Q}}(E_\tau^\times)$
  (the zero case is solved by $(y,r)=(0,0)$ and is not absorbable); cube
  pullback $c=h(a,b,z^3)$ with $v_w(c)=2v_w(a)+6v_w(z)-2\ge4$.
- **CONDITIONAL (no record claimed):** L6 witness-tie $a=1+2s$ gives a
  universal definition of $\mathbb{Z}$ in $\mathbb{Q}$ with **6** quantifiers
  *iff* the explicitly scoped global norm-assembly lemma (THEOREMS.md §L6)
  holds. Unconditional records stand: $\forall_{10}$ refereed (Daans
  arXiv:2301.02107), $\forall_7$ unrefereed (Sun arXiv:2607.28606).
- **EVIDENCE (bounded, not proof):** 61/61 guarded canonical global tied
  certificates, $A$ a nonsquare $w$-unit, $\Delta=\{2,w\}$ exactly,
  32 $\tau{=}0$ / 29 $\tau{=}2a/A$. An earlier unguarded table was discarded
  (it accepted nonunit $A=5$ at $w=5$); those buckets sit in `data/superseded/`.
- **L7 STRUCTURE (corollaries + bounded probe, 2026-08-15):** at odd places
  with $v(A)=v(B)=0$ the tie's exact cost is the wild places of $M_\tau$
  (L7a); coefficient-bad places ($v(A)\ne0$ or $v(B)\ne0$) carry real
  machine-witnessed obstructions (incl. a cancellation case with
  $v_5(AB)=0$), so the support is
  $\{2,\infty\}\cup\mathrm{wild}(M_\tau)\cup\{v:v(A)\ne0\text{ or }v(B)\ne0\}$; no
  fixed $(s,b,\tau)$ has a rational-function witness $Y,R\in\mathbb{Q}(z)$
  (L7b, conditional on L6 soundness), and probed canonical pairs certifiably
  miss parts of $\mathfrak m_w$; obstruction sets are always even (L7c).
  Probe: 93/93 target cells completeness-positive; non-target: 0 solutions
  found on 224 cells (205 certified failures, 19 refusals). Assembly stays
  open.
- **L8 THE ALIGNMENT WALL (2026-08-15):** exact value identities
  $(Ab^2D_z)^2u_\tau=X^2-D_\tau Y^2$ with
  $D_\tau\operatorname{disc}(E_\tau)\equiv-2AbP$ in both branches.
  **Proved:** absorption (odd $v$ with $v(A),v(B),v(M_\tau)$ all even is
  soluble — support shrinks into the odd-valuation places); the 2-adic
  parity wall $v_2(-2AbP)=1$ on every admissible witness (the source and
  target square classes $D_\tau$, $\operatorname{disc}E_\tau$ never
  coincide — the source algebra can be split when $P\in\mathbb Q^{\times2}$,
  while $E_\tau$ is always a field); emptiness of both exact
  square-class matching gauges (the first Pell-parametrized for $s\ne0$
  and swept; the second, and the $s=0$ slice, by 2-adic parity). **Evidence
  (scoped):** all 10 frozen rescues are wild-prime alignment events
  (largest odd-mult primes $1.6\times10^5$–$4.5\times10^{21}$, symbols all
  $+1$); recorded sweep 13/161 pairs soluble ($\approx8\%$, grid in
  NOTES). Assembly stays open.
- **L9 RECIPROCITY STEERING (2026-08-16):** **Proved:** the steering lemma
  (the product of $(x,d)_v$ over $T=\{2,\infty\}\cup\operatorname{supp}x
  \cup\operatorname{supp}d$ is $1$, so one unchecked place is forced) and
  the prime-$b$ symbol law ($b=\varepsilon q_1$ contributes exactly
  $(x,d)_{q_1}=(A|q_1)$ with $v_{q_1}(x)=-4$, under unit hypotheses
  imposed factor by factor — a product-form guard is false and is
  re-refuted in-suite; admissibility forces
  $A\equiv5\bmod8$, never square). **Evidence (bounded):** steered search
  covered **345/353 cells** (odd prime $w<100$; fixed 15-value $u$-pool
  filtered by $v_w(u)\ge0$) with proved-soluble
  admissible witnesses — 344 with the wild symbol *concluded* from the
  other places by the lemma and cross-replayed (largest forced prime 45
  digits); frozen as `_L9_STEERED`, full table re-verified every run. The
  8 residual cells are named search-exhaustions with spread bad-place
  profiles (no local wall). Fixed-pair baseline: 13/161. Assembly stays
  open; the count 6 stays conditional.

- **L11 BRANCH COMPLETION + THE BRANCH-INDEPENDENT WALL (2026-08-16):**
  **Proved:** adjoining a branch costs nothing (soundness is $\tau$-uniform; $k$
  conics in the same $(y,r)$ multiply to one polynomial) given a denominator
  nonvanishing on $\Phi$; the four square classes of $\langle-1,A\rangle$ are hit
  by explicit rational $\tau$; the controlled set
  $\{2,3,5,7\}\cup\operatorname{supp}(\alpha)\cup\operatorname{supp}(\delta)$ is
  complete. **The decisive result is negative:** the pairing
  $(x,d)_A(x,d)_q=-1$ holds on *every* branch when $v_A(z)\le0$, so branch choice
  cannot rescue a cell. Reachable **iff $z$ has a numerator prime $\equiv1\bmod4$**:
  **190 of 353 cells**, constructively certified; the other **163 are
  unreachable for the prime-$b$ family on every branch** — proved for prime $A$,
  evidence-only (878,400 certificates, 0 hits) for composite or rational $A$. An intermediate 353/353 claim was
  **false and is retracted**, and the same defect invalidated 130 of L10c's 164
  rows — both found by adversarial audit, both frozen as regressions. Count 6
  and efd 5 unchanged; assembly still open.
- **L12 THE GENERALIZED WALL (2026-08-16):** the obvious repair of L11g — a
  cofactor $m$ with $(m\mid A)=-1$ — **fails**, for a reason that generalizes.
  At every odd place of $b$ coprime to the cell data the unit part of $x$ is
  $\equiv A$, so $(x,d)_p=(p\mid A)$ and
  $\prod_{\{A\}\cup\operatorname{oddsupp}(b)}(x,d)_p=-1$ **for every admissible
  $b$** — prime, semiprime, squarefull, rational. The $-1$ relocates, never
  cancels; **the whole coprime-$b$ direction is closed on every branch**
  (231 in-suite instances; subsumes L11g). The single surviving hypothesis is
  coprimality: $b=\varepsilon fq$ with $f\mid z$ lights up **103 of the 163
  walled cells** — class alignment now **VERIFIED** (L13a: $S_{\min}=
  \{2,3,5,7,f\}$, 883 members, 0 violations).
- **L14 H ON THE GRID COMPLETE (2026-08-18):** hypothesis H is verified
  instance-wise on the entire 353-cell grid — **every one of the 293
  aligned classes (190 L11 + 103 ESC) carries a certified emergent-free
  member** (`_L13_H_CLASSES`, replayed every run: default 60-seed,
  extended all 293, both suites exit 0). Method: L13f cofactor ladder +
  BLS-relaxed Pocklington engine upgrade + three scan waves (frozen
  classes, upgraded-prover rerun, alternate aligned classes). The
  uniform all-$w$ statement of H remains the open analytic frontier.
- **L13 AFTERMATH (2026-08-17):** the class side is verified and the
  shortcut side is dead. (a) All 103 escape rows carry verified class
  alignment (minimal controlled set $\{2,3,5,7,f\}$, median modulus
  $2.4\times10^7$, 883 prime members, zero violations; emergent places stay
  step (ii) with exact even parity). (b) Grid **witness coverage 353/353 —
  COMPLETE**: 60 walled-no-escape cells fully soluble ($b=\varepsilon fq$
  sharing a prime with the cell data — $f\mid z$ in 57 rows, $q=w$ in 3;
  55 Hasse-Minkowski + 5 steered certs) plus 7 frozen residual-cell
  witnesses; the last cell $(89,(-2,1))$ closed 2026-08-18 by **L13f**:
  the cofactor-ladder zero-bad test (remainder of $P(b)$ after stripping
  $S\cup\operatorname{supp}(b)$ is a square OR a single proved prime —
  parity law forces its symbol $+1$; no full factorization) found 325
  certified soluble candidates across two boxes, 236 of them in the
  structured box; the frozen witness $(a,b)=(3,\,89/367)$, $\tau=19/37$,
  was lead-replayed rung by rung and is suite-asserted. (c) **No shortcut
  exists**: $P$ irreducible degree 8 over $\mathbb Q$ on all 706
  (cell, branch) rows (sympy + in-suite Frobenius certificates); content
  sqf $\in\{\pm1,\pm5\}$; 137.5M square-class tests and 40 designed
  families found **zero** hits — step (ii) cannot be collapsed
  structurally. (d) Emergent-free members exist: 5 certified zero-bad
  members, 9 found by sampling among 44 factored. Assembly hypothesis H
  (an emergent-free member of each aligned CLASS, uniformly) stays open;
  the count 6 stays conditional.
- **L10 CLASS-SIDE DICHOTOMY (2026-08-16):** the reduction's class step is
  settled per branch. **Proved:** the modulus freezes only the fixed places
  $S_0=\{2\}\cup\operatorname{supp}(\alpha)$, and the symbols that must all
  be $+1$ are $\Sigma=\{\infty\}\cup S_0\cup\{q_1\}$ (the moving prime $q_1$
  is held by L9a and its residue class, never by $N$; at every other place
  $v_p(d)=0$, so an even $v_p(x)$ gives symbol $+1$) — so the big fixed
  primes of $D_z$ are wild candidates, not
  controlled places; on the canonical branch $\tau=2a/A$ one has
  $\delta_\tau=1/A$ and $\alpha=-1$ **exactly**; and 164 aligned classes
  (exactly the canonical $w\equiv1\bmod4$ cells) are exhibited with $\gcd(q_1,N)=1$, so
  Dirichlet populates them — step (i) of the reduction is **done** there.
  **Wall (proved for prime $A$):** on $\tau=0$,
  $(x,d)_A\cdot(A|q_1)=(-2\varepsilon|A)=-1$ for both signs since
  $A\equiv5\bmod8$ (which requires $a$ odd — admissibility is load-bearing,
  $a=2$ breaks it), so no aligned class exists; evidence only for
  rational/non-squarefree $A$, and nothing is claimed uniformly in $s$.
  **Superseded by L11g/L11h:** the "164 of 353" count is wrong (130 rows were
  invalid) and the $w\bmod4$ split is not the right invariant — the reachable
  cells are the 190 whose $z$ has a numerator prime $\equiv1\bmod4$. L10b
  itself stands, as the $\tau=0$ instance of L11g.
  Assembly stays open; the count 6 stays conditional.

## Verify

```sh
cd math/h10q
python3 h10q.py                                        # full suite, ~25 s, must exit 0
python3 h10q.py --extended                             # + 61-prime tables, 6494 L5 pairs, L7/L9 full scopes, ~90 s
python3 l6_search.py --canonical /tmp/smoke.jsonl 3 5  # guarded searcher smoke
```

Dependency-free stdlib Python. All arithmetic runs on a proven-primality
engine: deterministic Miller–Rabin below the exact A014233 13-base bound,
generalized Pocklington certificates above it, explicit
`PrimalityBound`/`FactorBudget` refusal otherwise — no probabilistic
acceptance anywhere.


## Reproducibility scope (2026-08-19)

Exactly these L17 artifacts are **repo-regenerable** end-to-end from
checked-in scripts: `l17_ratemodel_padic.py` → `l17_badroots_closures.py`
(`data/l17_badroots_closures.jsonl`) → `l17_matched.py`
(`data/l17_matched.jsonl`); `l17_classrisk.py`
(`data/l17_classfactors.jsonl`); `l17_schinzel_audit.py`
(`data/l17_schinzel_audit.jsonl`); `l17_horizon{101,103,107,109}.py`
(horizon closures) and `l17_horizon109_crossk.py` (cross-k scan);
probes `l16_emergent.py`, `l14_replay_all.py`. The older L17 JSONLs
(`l17_sieve.jsonl`, `l17_badroots.jsonl` cohort, `l17_ratemodel_censored.jsonl`,
`l17_cofactor.jsonl`, `l17_cofactor_nonclosure.jsonl`) are **persisted
evidence without checked-in producers** — their data was reviewed and is
downstream-consumed, but regenerating them requires the original
session-side scripts (not recovered into the repo).

## Next (updated 2026-08-19)

The grid program is **closed**: all 353/353 cells witnessed; all 293/293
classes carry a verified emergent-free member (L14); the member-level
obstruction is matched-modeled exactly (L17, marginal 0.9993, union
$z=-0.09$); the Schinzel reduction's local conditions are audited clean
on all 293 canonical pairs (Dirichlet side proved, degree-8 prime-values
side open); the verified horizon now extends off-grid to
$w\in\{101,103,107,109\}$.

Frontier items:
1. **Horizon wave 2**: close $w\in\{113,127,131,137,139,149\}$ on the
   established off-grid machinery (`l17_horizon*.py` recipe); measure
   off-grid density against the validated model.
2. **Composite-$w$ frontier**: extend the machinery beyond prime $w$
   (L13e already covers the wall for composite squarefree $A$).
3. **The uniform step (ii)**: a density/independence statement over the
   prime-$Q$ subsequence (the emergent-freeness second factor of H);
   the rate band $c\in[0.017,0.80]$ per family brackets it with the
   refusal mass as the named dominant uncertainty.
4. Watch arXiv:2607.28606 for v2 (its §3 dyadic freezing and §7 norm
   approximation are the referee items; L5 here offers the exact
   identity for a corrected §5).
