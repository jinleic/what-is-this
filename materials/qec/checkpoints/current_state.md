# Checkpoint — 2026-08-22 (EXP-056 exact $[[162,8,14]]$; fixed-point screen closed through $n=162$)

**Exact distance method.** The reciprocal-pole transversal now feeds a complete
class-orbit SAT certificate. On the Wang–Mueller $(3,27)$ constructor
$A=1+y^{10}+y^{14}$, $B=y^{12}+x+x^2$, the $255$ nonzero logical classes
partition into **20** orbits under the verified order-$162$
translation/$x$-reflection group (sizes $6,9,18$). Each representative is an
affine coset of 77 independent sparse $Z$-stabilizer rows. Every $H_X$ column
has odd degree three, so every kernel word is even and cap 12 excludes through
13.

**Exact result.** Kissat returns UNSAT on all 20 class CNFs and repeats all 20
on fresh digest/version-bound replay. Initial serial solver time is 1,312.75 s
(39.79–89.19 s, median 62.00 s); replay is 1,287.67 s. An explicit weight-14
word passes NumPy and bitset kernel/non-stabilizer checks on both sectors under
the exact BB duality permutation. Therefore
**$[[162,8,14]]$ with $d_X=d_Z=14$ exactly**. Canonical certificate:
`results/certificates/exp056_wm_162_8_14_distance.json`.

**Soundness repairs.** Generic fixed-functional sectors must drop the
monolithic origin anchor (FR-027; a two-variable counterexample changes SAT to
UNSAT). The canonical result uses only the class route with verified
class-stabilizer anchors; the two-generator pole route is diagnostic and
cannot certify. Any verified SAT or non-replay-complete assembly atomically
revokes the authoritative certificate before writing convenience evidence.
EXP-055 admission invokes a pure EXP-056 validator that rebuilds the current
constructor, matrices, orbit cover, 20 CNF hashes/replays, duality and upper
witness from the same inode-bound snapshot it fingerprints. Hostile rereview:
**NO_BLOCKER**.

**Literature audit.** The reciprocal-pole isomorphism remains 27/27; 25 rows
reproduce printed $k$. Six literature rows are now local two-sided exact
certificates, zero ceiling violations, slack min/median/max 2/20/22. The source
$[[162,8,14]]$ value remains labelled a BP-OSD estimate; only the EXP-056
certificate promotes it to a threshold.

**Fixed-point screen, complete through $n=162$.** Thirteen nonempty-frontier
odd lattices; weight-3 $A,B$; $8\le k\le24$; 2,132 symmetry classes
representing 51,769 normalised pairs. All **1,928/1,928** referenced classes
are dominated: 1,804 reduced-pole witnesses, 119 CDCL witnesses, five exact
CP-SAT fallbacks. All 1,923 persisted witnesses are rechecked physically. 204
high-$k$ classes have no reference; zero survivors, zero undecided. Shards bind
the census, certificate-hashed references, pure-validator version and solver
protocol.

**Scope.** Algebraic $k$ census: all 65 odd lattices through $n=360$
(4,229,823,962 pairs). Exact fixed-point distance screen: through $n=162$.
The next actionable extension is to hash-bind the exact $n=180$ CSS frontier
before screening the nonempty frontiers at $n=170,186,198$.

Verification: 49 targeted odd-lattice/claims guards; full suite **972 passed, 1
skipped (973 collected)** in 796.48 s; hostile review **NO_BLOCKER**; paper
rebuilt to **17 pages / 708,817 bytes** with zero undefined references and no
oversized float; bundle `dist/pbb_nogo_bundle_2026-08-22.zip` verified
(1,704 zip entries, 0 SHA-256 mismatches).

# Checkpoint — 2026-08-21b (the demote law is an ideal invariant; trinomials can never be mixed; J-A resolved)

**Theorems J-G/J-H/J-I/J-J (proved; `notes/theorem_jg_ideal_invariant.md`, ledger 51–52).**
The single-row demote law is a function of the ideal $I=\operatorname{Ann}_R(A,B)$ alone:
$\dim S=2\dim I^\infty$, with demote-full $\iff I$ nilpotent, immune $\iff I^2=I$,
mixed $\iff 0\neq I^\infty\neq I$. EXP-053 reproduces EXP-052's 202-parent
classification from the ideal in **8.9 s, zero mismatches**, and sharpens it:
$I^2=0$ on all 192 demoting parents, $I^2=I$ on all 10 immune. **Odd lattices
cannot demote at all** (semisimple $R$; 3,600 parents checked) — among the seven
published BB instances only $[[90,8,10]]$ on $(15,3)$ is protected this way.
Splitting $G=G_2\times G_{\rm odd}$ gives the exact-vanishing coset law: exact
vanishing costs two support points per occupied 2-coset, hence **no
weight-$\le3$ pair is ever mixed on any lattice** (the whole literature family),
while weight 4 attains it — $A=(1+x)(y+y^2)$, $B=Ay$ on $(2,3)$ is mixed with
240/255 classes demoting, verified by ideal chain, module chain and brute force.
**Problem J-A resolved.** EXP-054 census: 653,022,021 weight-$\le3$ pairs over 18
lattices (all catalogue lattices + $(12,12)$, $(15,12)$ + a parity grid), zero
mixed, 450 independent cross-checks. Coset criterion separates the catalogue
exactly (10 single-coset immune vs 192 multi-coset demoting).
Suite: 13 new checks + 3 claims guards; paper 14 pp.

# Checkpoint — 2026-08-20 (J-E + J-E′′ machine law; X-side fully classified; envelope through n=180 dominated)

**Theorem J-E.** Demotion-immunity of a parent ⟺ $1\in L_{\rm pre}+\operatorname{Ann}_R(M)$
($M=\ker H_X/S_Z$, $L_{\rm pre}=$ left-null $[A\;B]$). Proof: demote-classes are
$\{y:y\notin L_{\rm pre}y\}$ by $R$-equivariance; (3)⇒(4) is Nakayama; (4)⇒(1) direct.
Demotion audit = one GF(2) rank test. EXP-047 (158 s + resume, exhaustive quotient
classes): **192 demote / 10 immune, exact**; invariant flag reproduces immunity on
202/202. **J-E′′ proved + independently audited, machine law exact on ALL 202 parents (EXP-052, ledger 50)**: demote fixed set = $I^\infty M = e_I M$; fractions exactly 1 on 192 demote-full and exactly 0 on 10 immune, **zero mixed**; chains ≤ 3, $I^2M=0$ on every demoting parent; the six $k_P\in\{24,40,60\}$ parents are all demote-full by the algebraic F-chain; enumeration census (EXP-050: 2,693,100 classes over $k_P\le20$) agrees to the integer on all 196. Universal no-mixed = open problem J-A. **X-side finalized via EXP-049+051:
192 demotion-realized (99 exact strict-collapse, 98 full-k) / 8 CERTIFIED
X-monotone (exhaustive channel-exclusion) / 2 X-U with verified weight-8
channels ($=202$).** All 133 distance-bounded EXP-046 D=0 siblings were rebuilt with identity
gates, independently verified, and finished terminal-UNSAT after bounded descent
(128 at cap 5, 5 at cap 3): 128 exact distance-6 and 5 exact distance-4 children.
Corollary J-E′ (λ-partners of a fixed seed; ledger 42).
Envelope refresh (ledgers 43/46): **within the 368-row published catalogue,
every row through $n=180$ is dominated** ($254/254$ at $n\le144$, $64/64$ at
$n=180$). The last holdout `phase2_107` (PBB exact $[[180,20,6]]$) was matched
with equality by its own certified CSS parent $[[180,20,6]]$ (fingerprint
`15654a1cebe3e639`); former leads `15_6_0224/0229` were earlier dominated by
the same certified CSS $[[180,8,16]]$ (`855e8bf...`). The universal envelope
statement — over ALL codes, not just the catalogue — remains **open**; only the
$n=360$ catalogue rows are unresolved (28 candidates + 22 pending).
X partition finalized 2026-08-20 via EXP-049 + EXP-051: 192 demotion-realized /
**8 CERTIFIED X-monotone** — `9_6_0175` + `phase2_75/76/77/83/84/87/109`;
exhaustive one-shot CP-SAT channel exclusion (dual-kernel membership encoding,
weights below $d_P$, INFEASIBLE certificate; the earlier basis+pairs scan is
superseded) / **2 X-U with verified weight-8 channels** (`15_6_0256`,
`30_6_0289`: $d_P=10$; Theorem-H family-closed ($T=8=k_P$ resp. $T=16=k_P$:
no PBB over them exceeds $d_Z(P)$ with $k>0$) and each carries an
independently verified weight-8 multi-row light $X$-stabilizer: monotonicity
undecidable by exclusion).  Provisional matched MC pair (ledger 44): Gross $6.33\times10^{-3}$ vs PBB
$2.21\times10^{-2}$, disjoint CIs; decode-core/shot ~10×; one cell of 10, not a verdict.

**Fleet**: exp039 n180/360 sweep at 139/202 parents certified (63 remaining, all $n\ge180$); matched-MC accumulates (PBB sched#5020 stage 13+, 56k shots). Escape leads 15_6_0224/15_6_0229 ride the sweep's n=180 frontier.

# Checkpoint — 2026-08-19b (Theorem J-C: collapse is decidable; boundary instance named)

**Theorem J-C (EXP-046).** The J.5-collapse mechanism is **decidable in polynomial
time**: for $c\in\ker A$, $D=0$, $C=\operatorname{circ}(c)$, validity ($M=0$) and
partner-membership are automatic; demotion ⟺ one GF(2) rank test
($z_0\notin S_Z+\Delta_C$); $\dim\ker A\ge k_P/2$. Machine-complete over **all 202
catalogue parents** (deterministic, 228 s): **192 demotions (all single D=0 basis
vectors), collapse on 99/100 bounded open-window parents, 98/99 k-preserving**.
Boundary instance `9_6_0175` (d=8): row-overlap $\le2$ blocks multi-row demotion;
single-row channel negative over $\ker A$ + full 58-dim syzygy basis + all pairs.
Zero mismatches vs sampled census. Ledger 40; paper §6 updated (md/tex/PDF, 12 pp,
0 undefined, claims lock 9/9).

**Sharp open instance**: is $\mu(\mathrm{Syz})\cap\operatorname{im}\Delta=\emptyset$
on `9_6_0175` — an exact quotient computation. 3 unbounded parents
(`30_6_0289`, `15_6_0256`, `phase2_109`) await the exp039 sweep for d-bounds.

# Checkpoint — 2026-08-19 (J.5 refuted; EXP-045 B′ decided; two-sector picture complete)

**J.5 — REFUTED (machine-verified, independently re-derived).** The X-side
monotonicity $d_X(Q)\ge d_X(P)$ is false: syzygy-family ($M=0$) perturbations demote
weight-6 parent X-stabilizers to genuine X-logicals of Q at full $k$ ($\dim\bar\Delta=0$).
Four verified witnesses: 6_6_0099 ($\le6{<}8$), 9_6_0136 ($\le6{<}12$), 12_6_0193
($\le6{<}12$, the flagship parent), 15_6_0219 ($\le6{<}14$). Inspection: $\dim V=112$,
$\dim\mathrm{Syz}\!=\!78$ → density $\sim2^{-34}$: why uniform hunts returned 0
violations. Lemmas: L1 ($d_Z(Q)\ge d_Z(P)$ **always**, proved), L3 (single-row demotion
$\Leftrightarrow$ $M=0$), L4 (syzygy excess $=k_P$, 900/900), L5 (no Z-mirror).
Main re-verified all predicates on a second GF(2) path; ledger row
`j5::xdistance_collapse_refutation` (REFUTED classification). Caps in Theorems G/H/I
untouched (pure-Z survivors only). Hazard for the FT program: a PBB sibling of the
headline [[144,12,12]] is silently [[144,12,≤6]] at full $k$.

**EXP-045 — Conjecture B′ decided on the only residual parent** (`9a7638586033`).
Over the ENTIRE valid-perturbation family ($\dim V=112$; exhaustive weight-≤6
enumeration = 2,533,006,645 instances, exact $\binom{112}{k}$ counts; 400 dense-draw +
154,602-pivot cross-checks clean): $\dim\bar\Delta\in\{0,2\}<4=T$ always → the
Theorem-H absorption premise is never met → **$d_Q\le6$ for every valid $[C\,D]$** of
this parent. Residual witness has no exceptional direction.

**Census in flight**: `CensusFinisher` (n=180/360 prevalence of the collapse family).
EXP-044 n180 log exists; agent died pre-harvest, replacement running. Sweep: 4 nice-0
shards at stride 4 (past the first post-respawn parents); MC at stage ~23 on
PBB sched#4107.

# Checkpoint — 2026-08-18 (W4 [5,12] certified both circuits; X-sector proof repair; EXP-042 emptiness; package artifacts)

**W4 circuit distance (EXP-040):** `certified_d_DEM_mech >= 5` on BOTH Gross
and PBB `[[144,12,12]]` circuits (flattened undecomposed structural DEM); ub 12
(EXP-029 weight-12 witnesses). w=3 exhaustive; w=4 meet-in-the-middle
COMPLETED on both — not budget-capped (`results/processed/exp040_circuit_distance.json`,
`verdict.certified_d_dem_mech_both_ge_5`). Intervals now [5,12] both — the
matched-MC arm (resumed, `exp040-mc`, 4 workers, stage-persistent) carries the
separation decision.

**X-sector repair (Theorem J note):** the paper's Theorem G(ii) justification
sentence was wrong ("X-sector of the quotient is unchanged"). The pure-X
centralizer SHRINKS (e.g. `phase2_58` 40→12). Correct supports: rank identity
for G(ii); Lemma J.3 (`dim(Xcen/S_X)=dim(Zcen/S_Z)=k` for every stabilizer
code); J.6 (`Xcen(Q)=Xcen(P)∩ker[C;D]`, rank bound). J.0: `d_X(Q)=d_Z(P)`
(reversal symmetry). J.4: demotion space feeds allowed dressing. OPEN:
unconditional `d_X(Q)≥d_X(P)` (J.5/J.7) — the §6 paragraph states it as open;
no cap uses it. `notes/theorem_j_xsector.md`; `tests/test_sector_identity.py`
(2 tests); md+tex rebuilt (11 pages, no undefined refs).

**EXP-042 residual emptiness:** NO BB parent with `T < k_P/2` on any lattice
`l*m ≤ 24` — 4,621 exhaustive symmetry-quotiented trinomial parents, verdict
`SMALL_CLASS_EMPTY` (`results/processed/exp042_residual_probe.json`).
Conjecture-B' counterexample hunt now running for `l*m ≤ 56` (ResidualProbe
follow-up). X-side counter-hunt (J.5-collapse search over 368-row valid
perturbations) running (XSectorAnalyst follow-up → `exp044`).

**368-row dimension-identity lock (HostileReview LOW-2 closed):**
`tests/test_dimension_identity_368.py` — 737 parametrized checks recomputing
`k_P − k_Q = dim(Δ̄)` from the raw v9 matrix, parent and PBB rebuilt
independently (no copied proofs-table text). Full theorem lock: **783 passed,
1 skipped**.

**Figures + package:** `reports/fig_rate_distance.pdf`,
`reports/fig_trade_law.pdf`; arXiv kit complete on disk:
`reports/arxiv_package.md`, `reports/arxiv_metadata.json`,
`notes/arxiv_submission_snapshot.md`, `notes/pending_validations.md`.

**Evening increments (also in `math/PROGRESS.md`):**
EXP-042-EXT closes the residual class through **ℓ·m ≤ 40** (10,645 parents, 43
shapes, zero residual; 4×8/5×8 rank-full empty); EXP-044 clears the X side on
**all 116 n≤144 parents** (0 decrease witnesses, 6.47M valid perturbations
after a contained uint8-filter fix; core GF(2) confirmed via `gf2.linalg`).
Hostile-review round 2: NO_CRITICAL; 2 HIGH fixed (G(ii) false "rowspace
extends" sentence → projection argument, machine-re-verified; PDF tab
corruption) + 3 MEDIUM + 5 LOW all applied; full suite **904 passed / 1
skipped / 905 collected**; PDF synced (12 pp, 0 undefined, rendered values
pdftotext-verified).

**Ops note:** broker restart consumed 3 subagent results without their files
landing; two regenerated (one running). exp040 processes were SIGSTOPped by
the broker — resumed via SIGCONT; exp040-mc respawned at nice 0 (thread-count
cap, single-digit threads keeps the QEC share ≤ 14 threads = 50 % of 28).

# Checkpoint — 2026-08-17 (Theorem H parent-level no-go + flagship paper)

**EXP-039 / Theorem H (2026-08-17)** — the family-wide layer.  The dressing
space $\Delta=\{\lambda[C\,D]:\lambda[A\,B]=0\}$ is an $R$-*submodule* of the
$Z$-sector (Lemma 1: the coefficient set is an ideal of
$R=\mathbb F_2[x,y]/(x^\ell{-}1,y^m{-}1)$); hence absorbing one minimum-weight
parent $Z$-logical absorbs its whole translation orbit.  With
$T(P)=\dim (M(P)+S_Z)/S_Z$ ($M$ = orbit-span of all minimum-weight
$Z$-logicals), a *terminating UNSAT* enumeration pins $T$ exactly, and
$$d_Q>d_Z(P)\ \Rightarrow\ k_Q\le k_P-T(P)\quad\text{for \emph{every} }[C\,D]\text{ of }P.$$
Machine results: exact $T$ for **134/202** parents — exhaustively all $117$ at
$n\le144$ (the $68$ remaining are all $n\in\{180,360\}$; sweep running,
figures monotone lower bounds), **61 family-closed** ($T=k_P$) including the
**Gross code** $A=x^3+y+y^2$, $B=y^3+x+x^2$ ($T=12$) and **all eleven**
catalogue $[[144,12,12]]$ parents; $249/368$ rows capped a priori with no
per-row search.  Saturation: all **7/7** certified reversals sit at
$k_Q=k_P-T$ with slack exactly zero, hypotheses replay-certified (gate
non-vacuous 7/7).  Envelope: **all seven** reversals are CSS-dominated at
equal $n$ by the strongest of **162** certified-exact CSS candidates
(`exp036_envelope_check.json`, schema v2 — scope extended to EXP-027's pair,
dominator canonicalised to max $kd^2$; FR-023 records the 5-of-7 scope defect
found while writing the paper).  Flagship write-up: `reports/paper_pbb_nogo.md`
(claims locked by `tests/test_paper_claims.py`, 8 checks).    **Theorem I (forced saturation).**  $\Delta$ is the image of the left kernel
  $L=\{\lambda:\lambda[A\,B]=0\}$ with $k_P=2\dim L$ (elementary; verified on all
  202 parents), so $\dim\bar\Delta\le k_P/2$ always; with Theorem H, every parent
  with $T\ge k_P/2$ (133/134 certified; exception `9a7638586033`, $k_P{=}12$, $T{=}4$)
  satisfies $d_Q>d_Z(P)\Rightarrow k_Q=k_P-T$ **exactly**.  EXP-040 probe: 26,898
  unseen perturbations, zero upper-law violations, 496/496 strict increases on the
  law with $\bar M=\bar\Delta$ containment; `tests/test_pbb_theorems.py`.
Suite **159
passed / 1 skipped**; ledger 32 entries.

# Checkpoint — 2026-08-16 (Theorem F + EXP-036 SAT closure campaign)

**EXP-038 / Theorem G (2026-08-17)** — the structural layer.  For parent
$[A\,B]/[B^T\,A^T]$ and its PBB, the dressing space
$\Delta=\{\lambda[C\,D]:\lambda[A\,B]=0\}$ satisfies $k_Q=k_P-\dim\Delta$
(**proved**, verified 368/368).  The PBB's pure-$Z$ centralizer equals the
parent's, so any minimum-weight parent $Z$-logical missed by $\Delta$ certifies
$d_Q\le d_Z(P)$ with no PBB search; contrapositively
$d_Q>d_Z(P)\Rightarrow k_Q\le k_P-t$.  All 7 certified reversals lie in the
$k$-halving class at margin **exactly 0**, and the criterion is provably silent
on each (falsification gate, `tests/test_pbb_survival.py`).  Per-row cost fell
785 s → 0.3-4 s.  All five EXP-036 reversals now carry **exact** PBB distances
($d=6$, $d=6$, $d=6$, $d=10$, $d=10$), cross-agreeing across two independent code
paths.  EXP-037: 259/368 envelope-dominated, two live leads at $n=180$.
Suite 129 passed; ledger 30 entries.


## Status
**Proved (all $\delta$):** $k$ can only fall relative to the parent, and a
same-$(n,k)$ CSS shadow's *$Z$-distance* caps $d_{\rm PBB}$.

**Proved ($\delta=0$ only, 213/368 codes incl. all 14 $[[144,12,12]]$):** the
parent CSS BB code weakly dominates in $[[n,k,d]]$.  For the 155 codes with
$\delta>0$ the strict-provenance audit (EXP-027) settles the shape: parent
dominates **66**, **2 certified reversals** (`phase2_58`/`phase2_60`:
$\delta{=}4$ perturbation of a $[[72,8,4]]$ parent gives non-CSS $[[72,4,6]]$
— distance strictly **increased**, both exact, double-verified), **87
undecided** at budget (that artifact is immutable).  **EXP-036** (2026-08-16)
turns the Theorem-F CDCL method on those 87: **24 more dominations and 5 new
certified reversals** are decided and replay-stamped (`phase2_71`,
`phase2_72` at $\delta=6$; `phase2_88`, `9_6_0183` at $\delta=2$;
`12_6_0217` at $\delta=4$, $n=144$) — running totals **107 dominations / 7
reversals / 41 undecided** (all remaining at $n\in\{180,360\}$; budget-laddered
sharded waves running).  All reversal codes remain inside the
CSS envelope at equal $n$, so the envelope stands (machine-checked: **all seven** certified reversals — EXP-036's five plus EXP-027's two — have an explicit dominating certified-exact CSS code at equal $n$, `results/processed/exp036_envelope_check.json`); "perturbations never buy
distance" is refuted for $\delta>0$ at three lattice sizes.  **Track B**
(below) extends the domination evidence to the previously excluded
$(12,12)$ and $(15,12)$ lattices.

**Proved, basis-independent (EXP-023, Theorem C6):** every catalogue PBB
$[[144,12,12]]$ needs depth $\ge 8$ against the Gross code's 7, for any
one-ancilla schedule measuring **any generating set**: no stabiliser element of
weight $\le7$ has a nonzero $X$-part (CP-SAT INFEASIBLE, completeness
certified, all 14), and the pure-$Z$ subgroup has rank exactly 66 $<$ 132, so
every basis contains a weight-$\ge8$ generator.  Gross control: $w^*=6$ exact.
EXP-022's decision procedure is superseded.  *(Incidence-counting
cannot make this unconditional: ancilla-ancilla parity forwarding refutes the
**premise** "one direct data gate per check incidence" — a circuit measures
$Z_0Z_1$ and $Z_0Z_2$ with one gate at $q_0$ where $\deg(q_0)=2$.  It uses $T=3$
layers $\ge\deg=2$, so it does **not** violate $T\ge\max_j\deg(j)$; the
inequality remains **open** for general multi-ancilla circuits.  Counting also
reaches only $\ge7$ for 4 of the 14 members, including the benchmarked
`12_6_0193` (1008 incidences / 144 qubits = 7 exactly), so counting cannot
prove the family-wide $\ge8$ at all; that rests on max check weight.  Open Q4
in `proofs/pbb_structure.md` records both gaps.)*

**Measured, consistent but not individually decisive:** among five uniformly
sampled minimum-depth schedules per code, all five PBB point estimates exceeded
all five CSS ones (Mann-Whitney $p=0.0040$; schedule-mean ratio 3.55×, Welch
$p=0.014$) — but the predeclared Bonferroni pointwise test **failed**, and
neither code's schedule was optimised.  Decoding was **6–17×** slower on an
idle machine (p50 10.3×, p95 16.6×, p99 6.3×); our proposed cause (BP
convergence) was measured and **refuted**.

## Environment
* Current execution host: macOS Darwin 25.5 arm64, Apple M4 Pro, 14 CPU
  cores, 48 GiB.  EXP-013's isolated timing artifact was collected on an
  earlier 28-logical-CPU host (`ncpu=28` is persisted); latency is not
  transferred across hosts.
* Python 3.13.14 in `qec-codesign/.venv` (rebuilt from `pyproject.toml` after the
  `.venv` and `.git` directories were removed externally).  The current
  campaign suite passes 45/45 tests; the earlier clean-reproduction gate passed
  its then-current 19/19 tests with numpy drifted 2.4.6 → 2.5.2.
* numpy 2.5.2 · scipy 1.18.0 · stim 1.16.0 · pymatching 2.4.0 · sinter 1.16.0
  · ldpc 2.4.1 · galois 0.4.11 · ortools 9.15.6755 · python-sat 1.9.dev13
  · qldpc 0.3.3 · networkx 3.6.1 · sympy 1.14.0.
* Pinned artifacts (`third_party/manifest.yaml`): BivariateBicycleCodes
  `fa77e33`, qcode-discovery `87cdf0c`, qLDPC `f7c33eb`, quits `0f36b8c`.

## Established

### Baseline reproduction (Gates 2, 3, 7)
* All 7 Bravyi BB instances rebuilt independently; commutation, rank, $k$,
  check weight (6), qubit degree (6), connectivity verified by **two
  independent GF(2) code paths**.
* Certified exact distances, **two-sided**: every sector proved `INFEASIBLE`
  below the value *and* an explicit minimum-weight logical exhibited and
  re-verified through the independent GF(2) bitset path (weight, kernel
  membership, non-membership in the stabilizer row space, `is_logical`).
  All four persisted in `results/certificates/` with their witness supports:
  $[[72,12,6]]\,6$ · $[[90,8,10]]\,10$ · $[[108,8,10]]\,10$ ·
  $[[144,12,12]]\,d_X{=}d_Z{=}12$.
  (The first version of these certificates proved only $d\ge D$ and claimed
  exactness with no witness — see FR-009; regenerated, values unchanged.)
* **Gate 7 passed**: an independent path that *transcribed* the published
  schedule from `decoder_setup.py` (rather than deriving one) agrees with our
  constraint-solver circuit on every count — depth 7, 288 qubits, 864 CX/round,
  0 detector events in 10 000 noiseless shots.

### Theory (`proofs/pbb_structure.md`)
Lemma 0 (validity ⇔ $AC^{\mathsf T}+BD^{\mathsf T}$ symmetric — agrees with the
primary source) · Props 1–3 ($k_{\rm PBB}=k_{\rm BB}-\delta$, verified 368/368)
· **Theorem 1 + Cor 1** (δ=0 ⇒ parent CSS code dominates; **213/368**, all 14
$[[144,12,12]]$) · Theorem 2 (CSS shadow, all δ) · Prop 4 ($d_X=d_Z$ for BB,
375/375) · Lemma C1 (one-ancilla correctness criterion) · **Theorem C3**
(criterion *forces* depth 7 within the **translation-invariant class**;
TI $T{=}6$ `INFEASIBLE`; unrestricted 6-exclusion external, ASC) · Theorem C4 +
**Theorem C6**
(PBB $[[144,12,12]]$ needs depth ≥ 8 for **every generating set** — schedule-class
independent AND basis-independent, EXP-023; $\mathrm{rank}(V_7)=66<132$ all 14) ·
Proposition C7 (three TI-refuting weight-9 PBB members have exact unrestricted
depth 9: persisted depth-$w_{\max}$ witnesses, structural checks, and 0
detector/observable events in 4000 noiseless 12-round Stim shots each; EXP-033).

### Target certification (EXP-034/035; updated 2026-08-15)
* **EXP-034 [proved + exact, target-only — Theorem E]:** `12_6_0193` is
  **genuinely non-CSS**: not CSS after arbitrary row operations
  ($\dim S_X{+}\dim S_Z = 20{+}66 = 86 < 132$, two independent GF(2) paths),
  after any qubit permutation (rank invariance proof), after any local-H
  assignment (affine GF(2) system infeasible, $0{=}1$ certificate persisted),
  **and after any independent per-qubit local Clifford**: nineteen parity
  masks plus one one-hot parity row XOR to $0{=}1$ over GF(2) — a
  solver-free linear certificate, machine-rechecked at write time, with
  CP-SAT independently INFEASIBLE as corroboration.  Permutation conjugates
  of local Cliffords are local Cliffords and CSS-ness is
  permutation-invariant, so the whole row-op × permutation × local-Clifford
  group is refuted.  Row space **indecomposable** under the same group.
  Non-vacuousness control: the LC-inequivalent $[[5,1,3]]$ perfect code
  stays linear-feasible (its refutation genuinely needs one-hot
  integrality).  `results/processed/exp034_target_equivalence.json`,
  `results/certificates/exp034_target_full_lc_decision.json`.
* **EXP-035 [proved + exact — Theorem F]:** independent target-distance
  certification is **complete**: $d(\texttt{12\_6\_0193}) = 12$, two-sided.
  Lower bound by strata — exact MITM exclusion of all weight-$\le5$
  centralizer elements; complete weight-6 classification (exactly 72
  weight-6 zero-syndrome vectors, all in $S$, set-equal to the 72 stored
  pure-$Z$ checks; $720 = 72\times10$ collision splits as an internal
  identity); weights 7–11 excluded by **four UNSAT sector proofs** from the
  independent PySAT/CaDiCaL backend (equisatisfiable CNF; 502/674/422/576 s
  single-threaded) over the machine-checked translation-orbit reduction
  (72 verified translations; transported dual rank 24 on two GF(2) paths).
  Upper bound by the independently re-verified weight-12 witness.  The
  racing CP-SAT backend proved nothing (sector 0 UNKNOWN across days) —
  CDCL closed what CP could not.  Collector is fail-closed: frozen
  backend/status proof pairs, mandatory solution re-verification, and
  cross-backend contradiction fail-stops, all regression-tested.
  Canonical: `results/certificates/pbb_12_6_0193_distance.json`; the
  partial artifact remains as historical evidence, superseded.
* **EXP-033 binding hardened (FR-021):** the noiseless-check artifact now
  stores/tests the canonical-input SHA-256 and per-row canonical-JSON slot
  hashes; canonical routing additionally pins the exact three-instance set and
  depth $=$ certified_depth $=$ $w_{\max}$ $= 9$; evidence regenerated.
* **EXP-024 fail-closed (FR-020):** the normal production CLI now requires the
  attested v2 gate *before* any target/schedule/circuit work and hard-blocks
  until a real v2 driver exists; canonical benchmark artifact remains a
  zero-byte placeholder; structural verdict `NEGATIVE_STRUCTURAL` unchanged.

### Circuit and decoder (matched, $p=0.002$, 12 rounds)
| | CSS Gross | non-CSS PBB `12_6_0193` |
|---|---|---|
| data / ancilla qubits | 144 / 144 | 144 / 144 |
| max check weight | 6 | 8 |
| mixed checks | 0 | 72 |
| certified depth | 7 (TI class; unrestricted 6-exclusion external) | 8 (class-free exact: lower bound + valid TI witness) |
| 2-qubit gates / round | 864 | 1008 |
| gate types | CX, CZ | CX, CY, CZ |
| valid min-depth TI schedules (exhaustive in TI class) | 8496 | 9968 |
| schedule-averaged LER | $8.90\times10^{-3}$ | $3.16\times10^{-2}$ (**3.55×**) |
| decode p50 / p95 / p99 (isolated) | 118 / 579 / 2159 ms | 1220 / 9614 / 13629 ms |
| two-sector hooks ($n{=}72$ analogue) | 0 % | 58.3 % |

Welch $p=0.014$; Mann–Whitney $U=25/25$, $p=0.0040$. The stricter Bonferroni
pointwise test **fails** and is reported as failing. Schedule choice alone moves
LER by 2.8× *within* each code, so no best-schedule claim is made.

### Track B — held-out lattices $(\ell,m)=(12,12)$ and $(15,12)$ (EXP-018/019)
**Exhaustive CSS dimension envelope** over the corrected Bravyi canonical family
(repeated monomials excluded, FR-011): $(12,12)$ 609 961 pairs, best $k=128$;
$(15,12)$ 1 219 834 pairs, best $k=160$.  QEC-relevant classes ($4\le k\le 24$)
are unchanged by the correction.  Bravyi's $[[288,12,18]]$ is recovered at $k=12$.

**Structural held-out test (GATE 10):** 2096 non-CSS PBB codes built from exact
commutation null spaces (weight $\le 2$, enumerated exhaustively by CP-SAT).
Zero violations of Prop. 3 ($k_{\rm PBB}=k_{\rm BB}-\delta$, $\delta\ge0$),
Prop. 4 (parent $d_X=d_Z$), the Theorem-2 shadow identity $k(Q')=k(Q)$, or P2.

**Distance domination (EXP-019):** 400+400 children of 10+10 QEC-relevant
parents ($8\le k\le 24$), weight-$\le 2$ perturbations enumerated completely on
every parent (two zero-count parents confirmed by brute force).  **Every code is
weakly dominated by its parent CSS BB code**, certified by an explicit witness
re-verified through `is_logical` + symplectic weight.  Parent distances exact:
$d\in\{4,6\}$ on $(12,12)$; $d\in\{6,8\}$ on $(15,12)$.

**Surface-code geometry and mapped circuit baseline (EXP-020/030):** first-order
same-$(k,d)$ totals include data + syndrome ancillas.  For codes with verified
circuits: $[[72,12,6]]$ total 144 (5.9× SC overhead), $[[144,12,12]]$ Gross
total 288 (12.0×); for codes without circuits, $[[288,12,18]]$ total 576
(13.5×) and $[[288,16,12]]$ total 576 (8.0×).  EXP-030 then ran a rotated
surface-$d=12$, matched-$k=12$ circuit baseline at numeric $p=0.002$:
$Z$ memory 5/200k and $X$ memory 4/200k LER/round on 3444 qubits, versus
Gross's five-schedule mean on 288 qubits.  Point estimate: about 30× lower
surface LER but 12× more qubits; at only 4–5 surface failures the ratio interval
is broad (approximately 6–140×).  Noise locations and decoder are **not**
identical, so this is a mapped comparison, not same-noise Pareto proof.

**Track B verdict (within scope):** no PBB code on either held-out lattice
exceeds the CSS rate-distance envelope at its own $(n,k)$.  Scope is 20 parents
and weight-$\le 2$ perturbations — not lattice-wide — but the structural
ceiling $k_{\rm PBB}\le k_{\rm BB}$ holds universally, so even untested codes
face a CSS parent at least as good in $k$.

## Falsified or downgraded — 21 entries, all in `notes/failed_routes.md`
1. **FR-001** our "catalogue distance refuted" claim — wrong nontriviality detector.
2. **FR-002** edge-colouring schedules — 1402 noiseless detector firings.
3. **FR-003** over-broad shadow campaign — superseded by an algebraic argument.
4. **FR-004** decoder timings under shared load — quarantined, re-measured.
5. **FR-005** our "83 unschedulable codes" claim — general schedules exist.
6. **FR-006** an erratum we nearly filed against a paper that was correct.
7. **FR-007** EXP-007 did not hold the circuit fixed — downgraded, replaced.
8. **FR-008** our BP-convergence explanation for the decoder slowdown —
   measured and refuted (0 % convergence for *both* codes); the latency
   measurement stands, the causal claim is withdrawn.
9. **FR-009** our four "exact" distance certificates proved only $d\ge D$ and
   contained no witness; solver semantics fixed, all four regenerated with
   independently verified witnesses, values unchanged.
10. **FR-010** wrong centralizer constraint in the held-out domination test
    (`code.H` instead of `lambda_swap(code.H)`); unsound verdicts deleted,
    regenerated with witness verification and a non-CSS regression test.
11. **FR-011** "exhaustive" enumeration over a family containing non-trinomials;
    stage-1 artifacts regenerated, representatives stored as terms not indices.
12. **FR-012** partial runs could clobber canonical artifacts — three
    successive guard defects; final guard is coverage-aware and exercised.
13. **FR-013** a verdict field that reported counterexamples as POSITIVE —
    corrected, regression-tested.
14. **FR-014** certified TI refutations misfiled as "undecided" — correcting
    this exposed the schedule-class separation ($T_{\rm TI}>13$ vs
    $T_{\rm unrestricted}=9$, EXP-033).
15. **FR-015** an aggregate coverage sentence multiplied counts and implied
    untested lattices were covered — wording now generated from per-stratum
    counters, regression-tested.
16. **FR-016** Theorem C3 and EXP-016 silently dropped their
    translation-invariant schedule-class scope — all claim surfaces corrected;
    unrestricted Gross depth-6 exclusion is cited externally, not claimed.
17. **FR-017** EXP-012 reversed catalogue distance-bound provenance and called
    five parents weaker — EXP-027 certifies two reversals and leaves three
    undecided.
18. **FR-018** EXP-026's 1.16-million-decode Cartesian grid failed its
    feasibility gate before any production row completed; protocol v2 uses a
    finite BP+OSD/BP+LSD set, resumable chunks, per-shot deadlines, and
    isolation-gated timing.
19. **FR-019** EXP-025 completed 70 Gross stage-1 schedules but zero PBB
    schedules before its unbounded decoder work unit stalled; v1 is
    development-only and symmetric, held-out v2 was frozen before PBB results.
20. **FR-020** EXP-024's frozen-v2 execution gate was dead code on the normal
    CLI; production now fails closed before structural work, and the canonical
    benchmark remains unexecuted.
21. **FR-021** EXP-033's noiseless evidence was not content-bound to its slot
    maps; canonical-input SHA-256 and per-row slot hashes are now persisted,
    gated, and regression-tested; checks regenerated.

## Not done
Track E rare-event estimators (ordinary Monte Carlo suffices at the measured
$10^{-3}$–$10^{-1}$ rates), general hardware mapping, exact physical-location
circuit distance at $n=144$, flag/cat families beyond EXP-024's specified
unverified 4+4 construction, and the general multi-ancilla depth/space-time
tradeoff.  Weight-$\ge3$ perturbations on held-out lattices and a same-noise,
same-decoder surface comparison also remain open.  EXP-024/025/026 production
runs and the four EXP-035 distance sectors are
tracked in `checkpoints/next_actions.md`.
