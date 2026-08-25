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
| [`THEOREMS.md`](THEOREMS.md) | Precise frontier statements (L1–L31, A1–A3), full proofs of everything proved here, the unchanged chain showing classical Schinzel H implies the conditional six-quantifier record, the L23–L24 unconditional barriers, the L25–L27 coupling/tie frontiers, L28–L29's reciprocal trace/lift closure, L30's lower-degree quartic section plus branch/trace rigidity, and L31's exact reciprocal member plus cube/dispersion/AP1 boundary. |
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
| [`l19_tauzero.py`](l19_tauzero.py) + [`data/l19_tauzero.jsonl`](data/l19_tauzero.jsonl) | L19 variable-$a$, $\tau=0$, $f=w$ class-side break of the canonical 5-wall, with exact symbol audits and a separately labelled bounded scan. Replay from `math/h10q/`: `python3 l19_tauzero.py`. |
| [`l19_classexist.py`](l19_classexist.py) + [`data/l19_classexist.jsonl`](data/l19_classexist.jsonl) | L20 uniform class-existence constructor for every cell, including the 353 grid certificates and the independent fixed-$a=1$ collision diagnostic. Replay from `math/h10q/`: `python3 l19_classexist.py`. |
| [`l19_cell179.py`](l19_cell179.py) + [`data/l19_cell179.jsonl`](data/l19_cell179.jsonl) | L19 search ledger and exact kernel replay closing the last tried off-grid cell, $w=179$, with every refusal kept non-evidentiary. Replay from `math/h10q/`: `python3 l19_cell179.py`. |
| [`l20_admissible.py`](l20_admissible.py) + [`data/l20_admissible.jsonl`](data/l20_admissible.jsonl) | L20 historical Schinzel-pair audit: conditions (b)–(d) and positive lead were proved there; its then-open irreducibility entry is superseded by L22. Replay from `math/h10q/`: `python3 l20_admissible.py`. |
| [`l21_reducible_locus.py`](l21_reducible_locus.py) + [`data/l21_reducible_locus.jsonl`](data/l21_reducible_locus.jsonl) | L21a–c reducible-locus, constructed-branch escape, and corrected-palindromic identities, with exact applicability controls. Replay from `math/h10q/`: `python3 l21_reducible_locus.py`. |
| [`l21_irred_generic.py`](l21_irred_generic.py) + [`data/l21_irred_generic.jsonl`](data/l21_irred_generic.jsonl) | L21d generic and 353-fiber irreducibility certificates, quantitative-HIT dodge, first-$a$ audit, and the cross-branch reducibility hunt. Replay from `math/h10q/`: `python3 l21_irred_generic.py`. |
| [`l21_irred_direct.py`](l21_irred_direct.py) + [`data/l21_irred_direct.jsonl`](data/l21_irred_direct.jsonl) | L21e direct 353-cell no-go audits for mod $w$ and Newton polygons, plus the $a=1$ quartic trace criterion and auxiliary exact certificates. Replay from `math/h10q/`: `python3 l21_irred_direct.py`. |
| [`l22_elimination.py`](l22_elimination.py) + [`data/l22_elimination.jsonl`](data/l22_elimination.jsonl) | **L22a — PROVED, PRIMARY CLOSURE:** for every fixed $Z\in\mathbb Q^\times$, the canonical-branch polynomial with $\tau^\dagger=(1+2a^2)/A$, $\delta=-4a^4/A$ is irreducible in $\mathbb Q(a)[b]$.  Symbolic endpoint elimination; finite scans remain EVIDENCE. Replay from the workspace root: `nice -n 19 python3 math/h10q/l22_elimination.py`. |
| [`l22_reciprocal_cube.py`](l22_reciprocal_cube.py) + [`data/l22_reciprocal_cube.jsonl`](data/l22_reciprocal_cube.jsonl) | **L22b — PROVED, INDEPENDENT CLOSURE:** on the explicitly named $a=1,A=5,s=0,\tau=3/5=\tau^\dagger,\delta=-4/5$ specialization, reciprocal trace plus complete $2$-isogeny descent proves irreducibility for every rational $Z\ne0$. Replay: `python3 l22_reciprocal_cube.py`. |
| [`l22_infinity.py`](l22_infinity.py) + [`data/l22_infinity.jsonl`](data/l22_infinity.jsonl) | **PROVED local side lemma; NOT the closure:** four quadratic Laurent clusters permit degrees $2,4,6$ locally; L22a performs the global elimination. Replay: `python3 l22_infinity.py`. |
| [`l22_factor_tuple.py`](l22_factor_tuple.py) + [`data/l22_factor_tuple.jsonl`](data/l22_factor_tuple.jsonl) | **PROVED criterion/no-go; CONDITIONAL prime values:** factor tuples can weaken whole-polynomial irreducibility when an all-plus residue is separately known, but parity controls only the product of signs. This route relocates the obligation and is not the closure. Replay: `python3 l22_factor_tuple.py`. |
| [`l22_square_branch.py`](l22_square_branch.py) + [`data/l22_square_branch.jsonl`](data/l22_square_branch.jsonl) | **PROVED off-chain algebra/class theorem; NOT APPLICABLE to the six-count:** the cell-dependent $\lambda$ selected by HIT ranges through the infinite L11c family and would require an additional witness. Never cite this artifact as the record closure. Replay: `python3 l22_square_branch.py`. |
| [`l22_fiber_geometry.py`](l22_fiber_geometry.py) + [`data/l22_fiber_geometry.jsonl`](data/l22_fiber_geometry.jsonl) | **PROVED Noether reduction; OPEN auxiliary nonvanishing:** the geometric reduction is exact, but its remaining polynomial lemma is not proved. L22a makes that lemma unnecessary for the conditional chain. Replay: `python3 l22_fiber_geometry.py`. |
| [`l23_half_sieve.py`](l23_half_sieve.py) + [`data/l23_half_sieve.jsonl`](data/l23_half_sieve.jsonl) | **L23a–b — PROVED small-prime-clean theorem and exact barrier:** dimension-$1/2$ semilinear sieve with BV level gives $\gg X/(\log X)^{3/2}$ members clean below $X^{0.49}$; it does not control the two-large-bad-divisor sector and is not a global member theorem. Replay from workspace root: `nice -n 19 python3 math/h10q/l23_half_sieve.py`. |
| [`l23_fibration.py`](l23_fibration.py) + [`data/l23_fibration.jsonl`](data/l23_fibration.jsonl) | **L23c–d — PROVED Capell/rank/Brauer results; member still OPEN:** the selected refined protocol has $2\theta$ nonsquare, so the fixed-cell conic bundle has non-split rank $10$ and $\operatorname{Br}(X)/\operatorname{Br}(\mathbb Q)=\mathbb Z/2$; no checked low-rank theorem applies. Replay: `nice -n 19 python3 math/h10q/l23_fibration.py`. |
| [`l23_norm_section.py`](l23_norm_section.py) + [`data/l23_norm_section.jsonl`](data/l23_norm_section.jsonl) | **L23e — PROVED scoped norm-section no-go:** exact matching and every coefficient-linear section miss $\Phi$; no polynomial-in-$b$ direct norm section exists; residual named curves remain OPEN. Replay: `nice -n 19 python3 math/h10q/l23_norm_section.py`. |
| [`l23_rational_section.py`](l23_rational_section.py) + [`data/l23_rational_section.jsonl`](data/l23_rational_section.jsonl) | **L23e — PROVED generic-section and constant-pullback no-go/reduction:** no generic $b$-line point; the first admissible factor ansatz reduces to genus-$4$ curves, whose per-cell points remain OPEN. Replay: `nice -n 19 python3 math/h10q/l23_rational_section.py`. |
| [`l23_squareclass.py`](l23_squareclass.py) + [`data/l23_squareclass.jsonl`](data/l23_squareclass.jsonl) | **L23e — PROVED dyadic squareclass walls:** $P=\pm2b\,\square$ is $\Phi$-empty, $P=A\,\mathrm{Norm}$ forces the wrong dyadic symbol, the fixed-$P$ square locus is genus $1$, and the reciprocal Capell locus is empty. Replay: `nice -n 19 python3 math/h10q/l23_squareclass.py`. |
| [`l23_multivar.py`](l23_multivar.py) + [`data/l23_multivar.jsonl`](data/l23_multivar.jsonl) | **L23f — PROVED tested-family complexity:** the $a(0)=0$ diagonal has an exact forced $Q^4$ square and a degree-$18$ squarefree representative with factor degrees $2,16$; no universal minimum-$21$ claim remains. Replay: `nice -n 19 python3 math/h10q/l23_multivar.py`. |
| [`l23_absorption.py`](l23_absorption.py) + [`data/l23_absorption.jsonl`](data/l23_absorption.jsonl) | **L23/L24 bridge — PROVED diagonal identities and finite-cover reduction:** discovers the two zero-cost self-couplings and their formal tied-rank-$\le1$ count; L24 proves both images $\Phi$-empty, so the nominal five-count is vacuous. Replay: `nice -n 19 python3 math/h10q/l23_absorption.py`. |
| [`l24_diagonal_geometry.py`](l24_diagonal_geometry.py) + [`data/l24_diagonal_geometry.jsonl`](data/l24_diagonal_geometry.jsonl) | **L24 — PRIMARY DYADIC CLOSURE + geometry:** proves both orientations empty on $\Phi$, integral degree-$8$ cover, $V_4\wr C_2$ closure of order $32$, genus-$35$ actual-cell slice, and no base-rational section. Replay: `nice -n 19 python3 math/h10q/l24_diagonal_geometry.py`. |
| [`l24_diagonal_arithmetic.py`](l24_diagonal_arithmetic.py) + [`data/l24_diagonal_arithmetic.jsonl`](data/l24_diagonal_arithmetic.jsonl) | **L24 — exact arithmetic tower:** replays the two eliminants, finite-flat rank $8$, both dyadic obstructions, monodromy, genus-$3$ $c$-line reduction and genus-$53$ actual-$z$ pullback. Replay: `nice -n 19 python3 math/h10q/l24_diagonal_arithmetic.py`. |
| [`l24_diagonal_local.py`](l24_diagonal_local.py) + [`data/l24_diagonal_local.jsonl`](data/l24_diagonal_local.jsonl) | **L24 — PROVED odd-place classification:** orientation II is locally soluble on the standard $k=1<v_w(Z)$ stratum for every odd target prime; the route is killed at $2$, not at the target. Replay: `nice -n 19 python3 math/h10q/l24_diagonal_local.py`. |
| [`l24_diagonal_search.py`](l24_diagonal_search.py) + [`data/l24_diagonal_search.jsonl`](data/l24_diagonal_search.jsonl) | **L24 — bounded corroboration only:** $370$ cells, $2{,}013{,}141$ guarded bases and $4{,}026{,}282$ orientation attempts, zero hits.  Every-cell emptiness comes from the mod-$16$ proof, never this scan. Replay: `nice -n 19 python3 math/h10q/l24_diagonal_search.py`. |
| [`l25_scaled_coupling.py`](l25_scaled_coupling.py) + [`data/l25_scaled_coupling.jsonl`](data/l25_scaled_coupling.jsonl) | **L25 — PROVED local escape, global route OPEN:** exact scaled eliminants; $(2X,\rho)$ and $(2X,\rho/2)$ are $2$-adically admissible exactly on even/odd $s$, together covering all of $\Phi$; the smooth $w\ge7$ formula lies on a freely chosen unit-$b$ stratum and is **not** an aligned-target theorem; bounded scans/searches labelled evidence only. Replay: `nice -n 19 python3 math/h10q/l25_scaled_coupling.py`. |
| [`l26_reciprocal_tie.py`](l26_reciprocal_tie.py) + [`data/l26_reciprocal_tie.jsonl`](data/l26_reciprocal_tie.jsonl) | **L26 — PROVED reciprocal tie/local structure, uniform member theorem OPEN:** the even pullback $Z=z^2$ and tie $\eta=Z(b+1)/(Db)$ make $P_{\rm rec}=b^4T(b+b^{-1})$ with quartic $T$; the tied dyadic symbol is automatically $+1$, standard W1 targets lift, and $(2b|p)=(2(u+2)|p)$.  L28 proves the trace field and bad-sign extension; L29 classifies the distinct lift; L31 proves one exact member at \(w=13\), not a per-target theorem. Replay: `nice -n 19 python3 math/h10q/l26_reciprocal_tie.py`. |
| [`l27_triangular_shear.py`](l27_triangular_shear.py) + [`data/l27_triangular_shear.jsonl`](data/l27_triangular_shear.jsonl) | **L27 — PROVED complete local shear + three scoped section no-gos, global route OPEN:** the one-piece coupling $(y,r)=(2X+28\rho,sX+\rho)$ covers every $2$-adic $s$-parity; an exact character-sum theorem plus $43$ small-prime rows gives a smooth point on every standard aligned odd target; one guarded real sample works.  The linear-$B$ cancellation section, $\lambda=\pm1,\pm A$, and $y=0$ are excluded exactly.  The full degree-$8$ cover, member theorem and five-count remain OPEN. Replay: `nice -n 19 python3 math/h10q/l27_triangular_shear.py`. |
| [`l28_trace_field.py`](l28_trace_field.py) + [`data/l28_trace_field.jsonl`](data/l28_trace_field.jsonl) | **L28 — PROVED uniform trace-field closure, uniform member theorem OPEN:** on every canonical L26 even-pullback parameter, the trace quartic is irreducible over $\mathbb Q_2$; the norm of $2(\theta+2)$ has squareclass $A\equiv5\bmod8$, so the bad-sign extension is always nontrivial and $T(v^2/2-2)$ is irreducible.  L29 supersedes the formerly open distinct-lift question; L31 supplies one exact member only. Replay: `nice -n 19 python3 math/h10q/l28_trace_field.py`. |
| [`l29_reciprocal_frontier.py`](l29_reciprocal_frontier.py) + [`data/l29_reciprocal_frontier.jsonl`](data/l29_reciprocal_frontier.jsonl) | **L29 — PROVED complete reciprocal-lift classification + exact reductions, uniform member theorem OPEN:** for square $Z$, $\theta^2-4$ is square exactly when $v_2(Z)=-2$ or $\ge2$.  Proof channels: strong Hensel, exhaustive Eisenstein-factor congruences, and an ordinary two-slope resultant Newton polygon.  The pullback $Z=8z^2/(1+z^2)$ preserves every positive odd valuation and forces the split local stratum.  The slice $b=w\rho^2$ freezes the norm field but has no generic section.  L31 promotes its \(w=13\) row to a proved member; no per-target theorem follows. Replay: `nice -n 19 python3 math/h10q/l29_reciprocal_frontier.py`. |
| [`l30_quartic_frontier.py`](l30_quartic_frontier.py) + [`data/l30_quartic_frontier.jsonl`](data/l30_quartic_frontier.jsonl) | **L30 — PROVED constant-two quartic + all-target selected fibres + rigidity reductions, global point OPEN:** $(y,r)=(2,sX+\rho)$ yields an exact quadratic in $m=\lambda^2$, hence degree $4$ rather than L27's degree $8$.  A uniform dyadic Hensel proof and character theorem close every odd $w\ge5$; the exact fibre $(a,b,z)=(5,3,3)$ closes $w=3$ by strong Hensel, while $(1,3,3)$ is an empty fixed-$a$ control.  One real stratum works.  General rational square-branch ties have no generic norm section, general $B$-independent shear cancellation is impossible for $s\ne0$, and the square-pullback trace-base bundle retains a separate reciprocal-lift conic.  The bounded zero-hit scan is EVIDENCE only. Replay: `nice -n 19 python3 math/h10q/l30_quartic_frontier.py`. |
| [`l31_frontier_push.py`](l31_frontier_push.py) + [`data/l31_frontier_push.jsonl`](data/l31_frontier_push.jsonl) | **L31 — one exact reciprocal member PROVED; other closures OPEN:** \(w=13,a=3,Z=169,q_0=1,b=13\) has a recursive Pocklington certificate and every Hilbert symbol \(+1\), with the mandatory lift \((\rho,\lambda)=(14/13,12)\).  The constant-\(c=-64/25\) section gives an exact rational \(H_2\)-point before the cube pullback and reduces cube compatibility to two genus-\(2\) curves.  The external Magma audit below proves both have only \(z=0\), closing this one section.  The analytic audit records the positive main term that signed cancellation cannot remove and proves AP1 equivalent to intermediate H. Replay: `nice -n 19 python3 math/h10q/l31_frontier_push.py`. |
| [`data/l31_magma_genus2.json`](data/l31_magma_genus2.json) | **L31 external-CAS completeness audit:** exact Magma V2.29-9 request, output, version, and intrinsic-help semantics.  `RationalPointsGenus2` returns only \((0:\pm8:1)\) with completeness `true` for both cube-compatibility curves, so the unique constant-\(c\) section has no guarded nonzero cube point.  Scoped external proof; no checked-in standalone Magma certificate and no claim about other L30 sections. |
| Engine note (2026-08-18) | Proven-primality engine strengthened: `_pocklington` now uses trial division to $10^6$ (early abort at $F\ge n^{1/3}$, exact integer cube root) + the **Brillhart–Lehmer–Selfridge relaxation** (CP Thm 4.1.5: $F\ge n^{1/3}$, two-factor discriminant test). Every verdict still exact; `PrimalityBound` refusals remain non-evidence. Latest suites: default $58.56$ s, extended $179.74$ s. |
| `data/l6_witnesses.jsonl`, `data/l9_steered.jsonl` | **Canonical serializations** of the authorities `_L6_WITNESSES` (61 rows) and `_L9_STEERED` (345 rows) in `h10q.py`. Byte-identical match required on every suite run — drift, absence, or any underived field fails the suite. Regenerate: `python3 h10q.py --export-evidence`. |
| `data/l9_steer_run.jsonl`, `data/l9_rescue_*.jsonl` | Raw steered-search run artifacts — search-side provenance only, never citable evidence. |
| `data/l13_filter_run2.json` | L13f stage-2 outcome (2026-08-18): full box spec, 22 counters with reason-coded refusals, 4,481 stored cofactor rows with per-row verdicts, 236 hit records, dense_box block. Authority for the cell-closure counters; the suite-frozen witness is the closing claim. |
| `l13_dense.py` + `data/l13_dense_smallb.jsonl` + `data/l13_dense_cofactorbig.jsonl` | DenseScan dense deep box: every odd $|b|<20001$, four branches, $a\in\{1,17,25,33\}$ — 319,968 candidates, 3,041 aligned (ALL $a{=}17$), 1,286 cofactor-big rows (handoff to the L13f ladder → 89 certified soluble). Exact counters only. |
| `data/litscout_h10q.md` | LitScout-2 (2026-08-18): source-verified citations for CONDITIONAL.md §4 — Krumm JTNB 28 (2016) pins (Thm 1.3, Props 3.4–3.5; degree 3 conditional on the elliptic Parity Conjecture, Prop 3.8), degree-8 ceiling beyond proved range on both sides, DDF arXiv:2102.06941v5 attribution confirmed, Crelle 495 (1998) author-list correction. Applied to CONDITIONAL.md same day. |
| `data/superseded/` | Inadmissible or superseded raw search buckets and the old rescue log — provenance only, never citable evidence. |

## State (2026-08-24)

- **PROVED here (machine-checked):** W0 canonical two-branch selector from the
  character identity $\sum_{a\in k}\chi(1+4a^2)=-1$ over every odd residue
  field including $\mathbb{F}_3$; W1 exact tied target-place criterion
  $\chi_w(-\delta_\tau A)=1$; W2 exact global criterion — the disjunction
  $M_\tau=0$ **or** $-\delta_\tau A M_\tau\in N_{E_\tau/\mathbb{Q}}(E_\tau^\times)$
  (the zero case is solved by $(y,r)=(0,0)$ and is not absorbable); cube
  pullback $c=h(a,b,z^3)$ with $v_w(c)=2v_w(a)+6v_w(z)-2\ge4$.
- **CONDITIONAL RECORD ESTABLISHED:** assuming classical Schinzel H
  alone, L19–L22 prove the intermediate per-cell H and complete the L6
  witness-tie architecture.  Hence $\mathbb Z$ has a definition over
  $\mathbb Q$ with **six universal quantifiers**, and
  $\operatorname{efd}_{\mathbb Q}(\mathbb Q\setminus\mathbb Z)\le5$
  (`THEOREMS.md`, L6 and L22d; `CONDITIONAL.md`, §§1–2).
  Classical Schinzel H is unproved; the unconditional anchors remain
  $\forall_{10}$ refereed (Daans arXiv:2301.02107) and $\forall_7$
  unrefereed (Sun arXiv:2607.28606), as sourced in `THEOREMS.md`, A1–A2.
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
- **L19–L20 CLASS EXISTENCE + SCHINZEL ADMISSIBILITY
  (historical precursor, completed by L22):** one aligned class exists
  for every cell, using $f=w$.  L20 proved conditions (b), (c), and (d)
  and the positive-leading part of (a); literal $S$-support of the
  content is false, while square-class support is the correct sufficient
  statement (`data/l19_classexist.jsonl`;
  `data/l20_admissible.jsonl`).  Its sole then-open item,
  irreducibility in condition (a), is now PROVED by L22.
- **L21 IRREDUCIBILITY PRECURSOR (SUPERSEDED AS FRONTIER):** blanket
  irreducibility is false off branch when $s=0$ and
  $\delta_\tau$ is square, but the canonical branch
  $\tau^\dagger=(1+2a^2)/A$ has
  $\delta_{\tau^\dagger}=-4a^4/A$.  L21 proved generic
  irreducibility and the recorded $353/353$ vertical fibres; L22
  upgrades this to every fixed nonzero rational $Z$
  (`data/l21_irred_generic.jsonl`; `/tmp/l22_elimination.md`).
- **L22 FIXED-FIBRE CLOSURE:** the symbolic factor elimination proves
  $P(a,Z,b)$ irreducible in $\mathbb Q(a)[b]$ for every fixed
  $Z\in\mathbb Q^\times$ on
  $\tau^\dagger=(1+2a^2)/A$, $\delta=-4a^4/A$
  (`/tmp/l22_elimination.md`; `data/l22_elimination.jsonl`).  The
  reciprocal specialization at
  $a=1,A=5,s=0,\tau=3/5,\delta=-4/5$ gives an independent all-$Z$
  proof (`/tmp/l22_reciprocal_cube.md`;
  `data/l22_reciprocal_cube.jsonl`).  Quantitative HIT in L20's odd
  character-admissible $a$-progression then selects a concrete
  irreducible $a$ (`THEOREMS.md`, L21d and L22).
- **L23 UNCONDITIONAL FRONTIER:** the selected linear/octic sequence has
  a genuine dimension-$1/2$ semilinear lower sieve: with
  $z=X^{49/100}$ it gives $\gg X/(\log X)^{3/2}$ prime members clean at
  every bad root prime below $z$.  Degree $8$ still permits up to $16$
  large factors, and reciprocity only forces an even bad count.  The
  first missing input is the parity-sensitive two-large-bad-divisor
  sector $p_1p_2>D_{\rm BV}$, equivalently $R_{\rm bad}\le1$.
- **L23 ALGEBRAIC BARRIERS:** the selected protocol cannot make
  $2\theta$ square in the octic field; its fixed-cell conic bundle has
  non-split rank $10$ and a $\mathbb Z/2$ Brauer class.  Exact norm
  matching, polynomial sections, elementary squareclass loci, fixed
  twists and the tested moving-$a$ families do not supply a member.
  Pre-existing even-$v_w(z)$ classes with $w\mid s$ remain OPEN
  uniformly, although the free L20/L22 refinement avoids that stratum.
- **L24 DIAGONAL CLOSURE:** both exact self-coupled covers are finite
  flat of degree $8$, but both images are uniformly empty on $\Phi$ over
  $\mathbb Q_2$.  Orientation I has coefficient valuations $(4,0,0)$;
  orientation II is killed by the complete mod-$16$ contradiction.
  Therefore the diagonal $\le5$ count is vacuous, not open.  This does
  not alter the fixed canonical six-count or remove Schinzel H.
- **L25 SCALED-COUPLING ESCAPE:** L24's unit scalings are not the whole
  self-coupled family.  The fixed type-I couplings $(y,r)=(2X,\rho)$ and
  $(2X,\rho/2)$ are $2$-adically admissible exactly on the even-$s$ and
  odd-$s$ strata of $\Phi$, respectively, and hence their disjunction
  has no uniform dyadic obstruction.  On a freely chosen unit-$b$
  unramified target stratum, an explicit smooth construction works for
  every odd $w\ge7$ (the streamlined proof does not lift at $w=3,5$, and
  aligned-target compatibility remains separate).  Two guarded bridge
  samples give the parity-wise real certificates.  No rational coupled
  point, member theorem, five-count, or unconditional record is claimed.
- **L26/L28 RECIPROCAL TRACE FIELD:** replacing the cube by the even
  pullback $Z=z^2$ and tying the third block coordinate to
  $\eta=Z(b+1)/(Db)$ gives
  $P_{\rm rec}(b)=b^4T(b+b^{-1})$ with $\deg T=4$.  On the canonical
  square branch the tied dyadic symbol is automatically $+1$, the
  standard guarded W1 target lifts, and the emergent character descends
  by $(2b|p)=(2(u+2)|p)$.  L28 now proves uniformly that $T$ is
  irreducible over $\mathbb Q_2$, that
  $N(2(\theta+2))\in A\mathbb Q_2^{\times2}$, and hence that the
  bad-sign extension is nontrivial.  L31 proves one globally good
  member at \(w=13\); a uniform/per-target member theorem remains OPEN.
- **L27 ONE-PIECE SHEAR:** the single coupling
  $(y,r)=(2X+28\rho,sX+\rho)$ is $2$-adically admissible for every
  $s$-parity.  A Weil character-sum argument, with all $43$ odd primes
  below $197$ exhausted exactly, gives a smooth point on the standard
  ramified L20 target stratum for every odd $w$; a guarded real sample
  also works.  Thus the L25 parity, $w=3,5$, and aligned-target local
  gaps are gone in this new shear.  A rational point on its
  bridge-specialized degree-$8$ cover remains OPEN, so no five-count or
  unconditional record follows.
- **L28 TRACE-FIELD CLOSURE:** after $w=u-2$, the trace quartic has a
  one-segment Newton polygon; its Ferrari resolvent excludes every
  quadratic-factor shape.  Since $A\equiv5\bmod8$, the norm squareclass
  proves $T(v^2/2-2)$ irreducible.  L28 closes the trace field and
  bad-character extension; L29 below closes the distinct lift.
- **L29 RECIPROCAL LIFT + COMPRESSION:** for even $t=v_2(Z)$,
  $\theta^2-4$ is square exactly for $t=-2$ or $t\ge2$.  The
  target-preserving map $Z=8z^2/(1+z^2)$ forces $t\ge2$ on every nonzero
  rational fibre.  The fixed-squareclass slice $b=w\rho^2$ reduces the
  member problem to one quadratic field but retains parity pairs.
- **L30 CONSTANT-TWO QUARTIC:** fixing
  \((y,r)=(2,sX+\rho)\) turns the coupled equation into an exact
  quadratic \(H_2(m)=0\) in \(m=\lambda^2\).  Its leading coefficient
  has dyadic valuation \(2\) on \(\Phi\), so this is genuinely an even
  quartic, not L27's octic.  Hensel solves every dyadic stratum; a
  two-quadratic character count gives fibre-regular rows at every
  \(w\ge5\); and \((a,b,z)=(5,3,3)\) gives an exact strong-Hensel
  \(w=3\) fibre.  The empty control \((1,3,3)\) proves fixed-\(a\)
  nonuniformity.  One real stratum works.  A rational global root and
  the remaining controlled places are OPEN.  L30 also proves generic
  branch/shear rigidity and separates the square-pullback trace-base
  auxiliary from its mandatory reciprocal lift.
- **L31 EXACT MEMBER + SHARPENED WALLS:** the mandatory fixed-field
  reciprocal slice now has one proved global member:
  \[
  (w,a,Z,q_0,b,\rho,\lambda)
  =(13,3,169,1,13,14/13,12).
  \]
  Its norm value is \(2^4Q/(3^2 37^3 83^2 1033^2)\), with
  \(Q=14082426920623718389\) prime by recursive Pocklington, and every
  Hilbert symbol against \(26\) is \(+1\).  This is one target, not a
  uniform member theorem.  On the quartic side,
  \(a=1,c=-64/25,\lambda=3,b=-3253/3125\) gives an exact rational root
  at \(Z=2033125/6101423\), but \(Z\) is not a cube and the shared prime
  \(3253\) leaves \(c\) a unit.  The required cube step reduces to two
  explicit genus-\(2\) curves; an external complete Magma calculation
  finds only \(z=0\) on both and closes this unique constant-\(c\)
  section.  Other quartic sections remain open.  Analytically, the unsigned two-large
  root-pair expansion has a positive main term; the correct target is a
  fixed-family Buchstab/Hilbert-detector asymptotic, not signed-error
  cancellation alone.  AP1 (\(R_{\rm bad}\le1\) in one selected class
  per cell) is equivalent to intermediate H by even parity, but remains
  OPEN.


**Current chain, exactly.** For every cell, L22 plus quantitative HIT
selects an irreducible admissible $a$; L20 supplies the aligned $f=w$
class and conditions (b)–(d), while L22 completes condition (a).  L19
supplies a member under classical Schinzel H.  L31 proves that AP1 is
equivalent to the intermediate member hypothesis but does not prove
AP1.  Therefore
`classical Schinzel H => intermediate H => Theorem C`
(`THEOREMS.md`, L22d and L31e; `CONDITIONAL.md`, §§1–2).  No per-cell
irreducibility certificate remains.

**Unconditional frontier, exactly.** L23's small-prime cleanliness is
insufficient for a uniform globally good member.  L26 supplies a
reciprocal six-count tie; L28–L29 close its trace field and local
reciprocal lift; L31 now proves one exact fixed-field member at \(w=13\)
but not a per-target theorem.  L30's degree-\(4\) equation has an exact
rational point on L31's unique constant-\(c\) bridge section before the
mandatory cube pullback; the cube-compatible controlled point remains
open.  The analytic target is a pointwise fixed-family
Buchstab/Hilbert-detector asymptotic whose first unavailable term is the
two-large-bad-divisor sector beyond BV.  Classical Schinzel H remains
the sole conjectural input.  AP1 and a target-specific cube-compatible
quartic point are the exact remaining closures.
- **L14 H ON THE GRID COMPLETE (2026-08-18):** hypothesis H is verified
  instance-wise on the entire 353-cell grid — **every one of the 293
  aligned classes (190 L11 + 103 ESC) carries a certified emergent-free
  member** (`_L13_H_CLASSES`, replayed every run: default 60-seed,
  extended all 293, both suites exit 0). Method: L13f cofactor ladder +
  BLS-relaxed Pocklington engine upgrade + three scan waves (frozen
  classes, upgraded-prover rerun, alternate aligned classes). This is the
  bounded-grid result; the L19–L22 implication above is the current frontier.
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
  members, 9 found by sampling among 44 factored. The stronger “member in
  every aligned class” statement is not required: L20 supplies one selected
  class per cell, L22 supplies irreducibility, and L19 gives member existence
  **CONDITIONAL** on classical Schinzel H (`THEOREMS.md`, L19–L22).
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
python3 h10q.py                                        # default suite; latest measured 58.56 s
python3 h10q.py --extended                             # full 61-prime/6494-pair scope; 179.74 s
python3 l26_reciprocal_tie.py                         # reciprocal identity/local replay
python3 l27_triangular_shear.py                        # dyadic + all-target local replay
python3 l28_trace_field.py                             # uniform trace/Capell irreducibility
python3 l29_reciprocal_frontier.py                      # lift classification + compressed pullback
python3 l30_quartic_frontier.py                         # quartic local theorem + rigidity
python3 l31_frontier_push.py                            # exact w=13 member + cube/dispersion/AP1 walls
```

The checked-in producers use dependency-free stdlib Python.  All their
arithmetic runs on a proven-primality engine: deterministic
Miller–Rabin below the exact A014233 13-base bound, generalized
Pocklington certificates above it, and explicit
`PrimalityBound`/`FactorBudget` refusal otherwise — no probabilistic
acceptance anywhere.  L31's scoped genus-\(2\) completeness result is
the one external exception: the exact Magma V2.29-9 request/output and
help semantics are recorded, but no standalone Magma certificate is
checked in.


## Reproducibility scope (2026-08-24)

All seven L23 producer/data pairs, all four L24 pairs, and the L25–L31
pairs in the inventory above are repo-regenerable.  Replay them from the
workspace root, **one process at a time**, with the displayed
`nice -n 19 python3 math/h10q/<producer>.py` commands.  They write,
respectively,
`data/l23_{half_sieve,fibration,norm_section,rational_section,squareclass,multivar,absorption}.jsonl`,
`data/l24_diagonal_{geometry,arithmetic,local,search}.jsonl`,
`data/l25_scaled_coupling.jsonl`, `data/l26_reciprocal_tie.jsonl`,
`data/l27_triangular_shear.jsonl`, `data/l28_trace_field.jsonl`,
`data/l29_reciprocal_frontier.jsonl`,
`data/l30_quartic_frontier.jsonl`, and
`data/l31_frontier_push.jsonl`.  Every producer is stdlib-only and
single-process; no pool is used.  `data/l31_magma_genus2.json` is the
explicit exception: a live external Magma audit with an exact replay
request, not a producer-regenerated artifact.  The remaining L29 fibres
and L30/L31 bounded searches retain their stated EVIDENCE scopes.

All six L22 artifacts are likewise repo-regenerable from their paired
scripts:
`l22_{elimination,reciprocal_cube,infinity,factor_tuple,square_branch,fiber_geometry}.py`
write
`data/l22_{elimination,reciprocal_cube,infinity,factor_tuple,square_branch,fiber_geometry}.jsonl`.
Their theorem labels and non-overclaim scopes are listed in the
inventory above; finite scans are never promoted to a uniform proof.

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

## Next (updated 2026-08-24)

The Schinzel-conditional chain remains closed on the fixed canonical
branch.  L31 proves one exact reciprocal member and closes its unique
constant-\(c\) quartic section after the cube; the uniform and full
quartic problems remain open.

1. **Primary reciprocal-tie target:** extend the proved \(w=13\)
   fixed-field member to every target.  The general mandatory slice is
   a sign-decorated degree-\(16\) binary-form norm problem; one exact
   member does not establish a density theorem.
2. **Primary coupled-cover target:** leave the externally closed
   \(c=-64/25\) section and construct a different target-specific
   rational root of \(H_2(\lambda^2)=0\) satisfying \(Z=z^3\) and every
   controlled place.
3. **Primary analytic target:** prove a pointwise fixed-family
   Buchstab/Hilbert-detector asymptotic with positive zero-bad constant.
   Its first unavailable term is the two-large sector beyond
   \(D_{\rm BV}\); signed character cancellation leaves a positive
   root-pair main term.
4. **Exact conjectural replacement:** prove AP1 for every cell.  Even
   Hilbert parity makes \(R_{\rm bad}\le1\) equivalent to intermediate
   H inside the L19–L22 protocol.
5. **Scope guard:** classical Schinzel H remains unproved and the sole
   conjectural input to the six-count; H10/$\mathbb Q$ is open.
6. **Optional side geometry:** the L22 Noether nonvanishing lemma and
   absolute rationality of the L24 total surfaces remain OPEN but are
   not chain obligations.
