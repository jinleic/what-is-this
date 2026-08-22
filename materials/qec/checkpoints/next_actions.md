# Next actions — ordered by expected information gain

## 2026-08-21b update — the X-sector boundary taxonomy is SOLVED; ranked plan moved to a note

The item "*a clean invariant separating the 10 immune parents from the 192*"
(2026-08-19b, point 3) is **closed**: the invariant is the ideal
$I=\operatorname{Ann}_R(A,B)$ — immune $\iff I^2=I$, demote-full $\iff I$
nilpotent (Theorem J-G), and combinatorially $\iff$ all $A$-monomials share one
$G_2$-coset and likewise for $B$ (Theorem J-I; for $m=6$: $y$-exponent parity).
Problem J-A ("does a mixed parent exist") is **resolved** affirmatively at
weight 4, and impossible at weight $\le3$ on every lattice.
Point 4 of the 2026-08-19c list (asymptotic demote-fraction estimate) is
**obsolete**: the fraction is exactly $1-\frac{2^{\dim S}-1}{2^{k_P}-1}$ with
$\dim S=2\dim I^\infty$, no estimate needed.

The ranked plan now lives in [`../notes/next_breakthroughs.md`](../notes/next_breakthroughs.md):

1. **A — odd-lattice PBB sweep** (theory-selected collapse-free region; startable
   with existing code, reuses EXP-054's invariant cache).
2. **B — multi-row channel theorem** (Fitting-ideal analogue of J-G; correlate
   against EXP-051's 8 INFEASIBLE + 2 witnesses first — hours, decisive).
3. **C — weight-4 mixed map** (gives a *usable* safe perturbation subspace $S$).
4. **E — fault-tolerant one-ancilla mixed-stabilizer circuit** (construct with
   flags/cat states, or prove the class-restricted no-go): the decisive
   end-to-end item.
5. **D — per-factor cyclic distance bounds** (unlocks the $n=360$ envelope rows;
   the naive idempotent-ideal bound is dead, see FR-024).
6. **F — rare-event estimator** to finish the matched-MC grid affordably.

Still-live operational items from the older lists: exp039 $n=180/360$ sweep
(63 parents uncertified), the $n=360$ envelope rows (28 candidates + 22
pending), and the matched-MC grid (1 of 10 cells).

## 2026-08-19c update — Theorem J-E landed; X residuals collapse to three bound-gated parents

`9_6_0175` is **closed** (provably X-monotone, Theorem J-E + overlap + light-scan);
the exact-characterization track is complete.
Fleet note 2026-08-19d: exp039 CLI changed — respawn convention is
`exp039_nogo_module.py run --ns 180 360 --stride 4 --offset {0..3}` (NOT `--ns ... --shard`). What remains:

1. **Fleet (exp039 sweep)**: certify the 66 uncertified parents; the X-U trio
   (`phase2_109`, `15_6_0256`, `30_6_0289`) resolves the second the bound lands
   (EXP-047 machinery in place: decision is one rank test at weight 8 vs d_P).
2. **Envelope**: the 26-record $n=180$ tail is dominated; remaining open surface:
   38 $n=180$ rows needing EXP-037 row records (run the row-record pass), 28 $n=360$
   candidates + 22 pending (sweep-dependent). Lead `15_6_0224/0229` episode CLOSED
   (both dominated by certified CSS $[[180,8,16]]`855e8bf...`; candidates 51/52
   retightened 18→14/16 → dominated).
3. **Matched-MC grid**: gross sched#6479 running (28k shots at last obs); 10 cells per
   p needed for the verdict assembly; pooled-latency comparison ONLY via the dedicated
   W5 latency arm (stage quantiles are not poolable).
4. Optional theorem polish: is demote-fraction → 1 generically (flagship = 4095/4095)?
   An asymptotic statement on $L_{\rm pre}$ orbit sizes would explain it
   ($(1-2^{-\dim L_{\rm pre}})^{k_P}$-style estimates).

## 2026-08-19b update — Theorem J-C landed; X-side residuals are now named instances

The X-collapse theory is closed *modulo named residuals*; the Z-side two-sector
program is complete (Theorem I + EXP-045). Ordered residuals:

1. **`9_6_0175` — the boundary instance.** Only bounded parent with an open
   window and no constructed collapse. Its row-overlap is ≤2, so any collapse
   must be single-row; that channel was probed through ker(A) + all 1,711
   syzygy basis/pairs of its 58-dim kernel with zero demotions. Sharp question:
   is $\mu(\mathrm{Syz})\cap\operatorname{im}\Delta$ empty — that is an exact
   quotient computation (not $2^{58}$), worth one focused pass.
2. **Certify d for the 3 unbounded non-collapsing parents** (`30_6_0289`,
   `15_6_0256`, `phase2_109`) — the exp039 sweep fleet owns this; once a bound
   exists the EXP-046 constructor decides collapse instantly (demotion is
   already verified=False, immune-to-depth; so the decision is one rank test).
3. **Boundary taxonomy (optional).** The 10 no-demote parents all have
   $\dim\ker A\in\{4,12\}$, but small-kernel is *not* explanatory: all 34
   nullity-2 parents demote. A clean invariant separating the 10 from the 192
   would give J-C a fully proved universality criterion.
4. **Residual-saturation mechanism** (carried from 2026-08-18/19): why do all
   dressing rows of `9a7638586033` live in a fixed 2-dim subspace (dim V→A map
   ≤2)? Parent-invariant formula wanted; pending with ResidualSaturation.
5. **End-to-end track** (owner: matched-MC + exp039 sweep): unchanged —
   Pareto verdict needs the MC completion + $n\ge180$ T-certs.

## 2026-08-17 update — Theorem H landed; program pivoted to closure + circuit

0. **[RUNNING] EXP-039 n=180/360 parent sweep.**  10 shards `exp039-f0..9`
   (nice -n 5).  134/202 parents exact — exhaustively all 117 at $n\le144$;
   remaining 68 are all $n\in\{180,360\}$ (one parent can take hours; slowest
   completed: 14,461 s).  Figures in `exp039_nogo_module.json` are monotone
   lower bounds; do not cite without re-running `assemble`.

1. **[RUNNING — sub-agents] 2026-08-17 parallel batch:**
   - **NOTE 2026-08-18: CircuitBuilder agent expired from registry (broker restart); its two supervised processes (exp040-mc, exp040-distance) persist under Main's ownership. Completed: W1 FINGERPRINTS_LOCKED (circuits + DEM sha256 exact), W2 all_passed (12 slot maps, TableauSimulator audit). W4-gross certified d_DEM_mech >= 5, ub 12. In flight: W4-PBB exclusion; W3 matched MC at 4 workers; W5/W6 to be run by Main on completion.**
   - `CircuitBuilder` (EXP-040): the decisive end-to-end experiment — verified
     detector semantics (TableauSimulator audit), familywise-pointwise LER at
     $p\in\{0.0015,0.002,0.003\}$, budget-capped $w{=}4$ circuit-distance
     exclusion, production-decoder arm.  Work order in
     `notes/` once landed (`notes/exp040_evidence_module.md`).
   - `SaturationProbe` (EXP-040-sat): Conjecture A — is $\dim\bar\Delta=T$
     forced for reversals?  Next open theory question in the program.
   - `TCrosscheck` (EXP-041): DONE — 12/12 in-budget parents match the SAT T
     exactly (incl. three $T{=}k_P$ closures); Gross over budget by design.
   - `HostileReview`, `FreshnessScout`, `ClaimsLedger`: referee audit,
     novelty check, claims matrix.
   - `PaperLaTeX`: DONE — `reports/paper_pbb_nogo.{tex,bib}` compiled clean
     (9 pp, 0 warnings, table-fidelity diff 0 mismatches).

2. **AFTER CircuitBuilder lands**: fold EXP-040 results into
   `reports/paper_pbb_nogo.md` §6 (the "not settled" paragraph) and decide
   whether the flagship claim upgrades to an end-to-end no-go statement.

3. **Conjecture A proof attempt** (if SaturationProbe returns SUPPORTED):
   saturation $\dim\bar\Delta=T(P)$ would turn Theorem H into an exact
   trade law.  Algebraic angle: $\Delta\supseteq\bar M$ forced by the
   survival criterion; equality would follow if $\Delta$ can never absorb a
   non-minimum-weight logical without also absorbing an extra orbit — test
   on small lattices first (EXP-040-sat sample).

## Immediate (queued, mechanical)
1. **EXP-038 / Theorem G — finish the structural sweep.**  The survival
   criterion (`src/qec_research/codes/pbb_survival.py`, proof in
   `proofs/pbb_structure.md` Thm G) decides domination by linear algebra once
   the parent's exact CSS distance is known, in 0.3-4 s per row.  It is gated on
   the same $n=180/360$ CSS frontier as EXP-037, so those two tracks share one
   bottleneck: **certify more $n=180$ and $n=360$ CSS BB parents** and both
   sweeps unblock together.
   a. `uv run python experiments/exp038_survival_criterion.py run --ns 108 144`
      — should close every remaining $n\le144$ row essentially for free.
   b. Shards `exp038-s{a..e}` cover EXP-036's 41 open rows; they fall through to
      the expensive climb where the parent is unpooled.
   c. Open theory question worth a proof attempt: is margin $t-\dim\Delta = 0$
      *forced* for every reversal, or a coincidence of these 7?  A proof would
      turn Theorem G (iv) into an exact characterisation of when perturbation can
      buy distance, and with it a no-go theorem for the $k_Q/k_P>1/2$ classes.
   Never claim a row dominated from the criterion alone when $d_Z(P)>d_X(P)$:
   the certified conclusion is $d_Q\le d_Z(P)$, not $d_Q\le d_P$.

2. **EXP-037 envelope classification — two live escape leads.**  The
   catalogue-wide question ("does any of the 368 published PBB codes escape the
   CSS BB envelope at its own length?") is 259/368 certified dominated with
   **two live leads**: `15_6_0224` (index 28) and `15_6_0229` (index 48), both
   $n=180$, $k=6$, verified witness weight 16 versus best certified same-length
   CSS $(k_C,d_C)=(8,14)$.  Decisive next steps, in order:
   a. `uv run python experiments/exp037_envelope_classification.py lower --indexes 28 48 --cap 14`
      — UNSAT at 14 proves $d\ge15>14$ and makes each a genuine escape
      candidate; SAT tightens $U\le14$ and dominates it immediately.
      (Dedicated probes `exp037-decide-28`/`-48` are running.)
   b. Complete the $n=180$ CSS frontier: 46 distinct parents, only a few
      certified so far.  A same-length CSS code with $k_C\ge6$ and $d_C\ge16$
      would close both leads without any PBB lower bound.
   c. Populate the $n=360$ CSS pool (39 distinct parents); the 14 flagged
      $n=360$ rows are pool-coverage artifacts, not findings, until it exists.
   Never report a flagged row as an escape without (a) a PBB lower bound and
   (b) frontier completeness at that length.

3. **DONE 2026-08-15 — EXP-035 exact distance certified: $d = 12$ (Theorem F).**
   All four orbit-representative sectors proved UNSAT by the independent
   PySAT/CaDiCaL backend (502/674/422/576 s single-threaded); canonical
   certificate at `results/certificates/pbb_12_6_0193_distance.json` with
   witness re-verification and transported dual rank 24 on two GF(2) paths.
   CP-SAT proved none of the four; the user approved stopping the redundant
   fleet and `pbb-s0..s3` are all stopped.
   The partial artifact is retained as superseded historical evidence.
4. **DONE 2026-08-15 — EXP-034 full-LC decided INEQUIVALENT.**  Primary
   proof is the 20-row GF(2) XOR contradiction (linear relaxation
   infeasible); CP-SAT independently INFEASIBLE.  See Theorem E in
   `proofs/pbb_structure.md` and
   `results/certificates/exp034_target_full_lc_decision.json`.
3. **Run exp013 on an idle machine** — the only admissible decoder-latency
   evidence.  It self-refuses if the 1-minute load average exceeds the
   threshold.  Command:
   `PYTHONPATH=src .venv/bin/python experiments/exp013_isolated_latency.py 400 12 0.002 4 4.0`
4. **EXP-024/025/026 v2 production** — all three are frozen and fail closed
   (FR-018/019/020).  Each needs its supervised, deadline-bounded, resumable
   driver implemented and its attested gate written before any claim-bearing
   run; EXP-024's normal CLI hard-blocks until then.

## High information gain
4. **Close the $\delta>0$ gap.**  exp012 on $n\le108$, $\delta>0$: 49 dominated,
   5 parent-weaker, 23 undecided (timeouts).  Re-run the 23 undecided with a
   longer budget; then decide the 5 "parent weaker" cases against the *certified*
   CSS reference set instead of the parent.  Target statement: a single number
   for "how many of the 368 published PBB codes are dominated in $[[n,k,d]]$ by
   an explicitly exhibited CSS code".
5. **Escapers.**  `phase2_58` and `phase2_60` ($[[72,4,6]]$, $\delta=4$) beat
   their own CSS shadow ($d_{\rm shadow}=4 < 6$).  They are still dominated by
   the certified Bravyi $[[72,12,6]]$ ($k$ 12 vs 4, same $d$) — but that needs
   to be stated as an explicit dominance claim, not left implicit.
6. **Y-gate PBB codes.**  Several catalogue codes have $A\cap C\ne\emptyset$, so
   their checks carry Pauli $Y$ and the max weight drops below
   $|A|+|B|+|C|+|D|$.  Two catalogue codes reach PBB weight 7.  Do any reach
   depth 7?  If one does, the depth separation for *that* code disappears and
   only gate count and hook structure remain — worth knowing precisely.

## Theory
7. **Even-depth pattern.**  Translation-invariant PBB schedules realise depths
   $8,10,12,14,16$ and never an odd value.  Derive (or refute) a parity
   invariant of Lemma C1 for group-structured codes.
8. **Hook-error asymmetry.**  For a CSS $X$-check an ancilla fault propagates
   only $X$; for a mixed check it propagates both $X$ and $Z$.  Formalise this
   as a lower bound on the circuit-distance loss of mixed checks — this is the
   mechanism behind the measured LER gap and is currently only argued verbally.

## Explicitly not attempted (state as such, do not imply otherwise)
* Verified flag-qubit or cat-state extraction beyond EXP-024's specified
  *unverified* two-ancilla 4+4 construction, whose structural verdict is
  `NEGATIVE_STRUCTURAL` and whose canonical Monte Carlo benchmark was never
  executed (FR-020; the depth bound in Theorem C4 is for **one-ancilla**
  circuits).
* Catalogue-wide local-Clifford equivalence classification.  EXP-034 settles
  the benchmark target **fully** (Theorem E: genuinely non-CSS under row
  operations, qubit permutations, local-H, and arbitrary independent local
  Cliffords, with row-space indecomposability under the same group); **any
  catalogue-wide statement remains open** — note the $[[5,1,3]]$ control
  shows the linear certificate does not decide every code.
* Rare-event / splitting estimators (Track E) — the measured logical error
  rates are in the $10^{-4}$–$10^{-2}$ range, where ordinary Monte Carlo has
  ample failures, so no rare-event machinery was required for these claims.
* New code discovery on the $(12,12)$ and $(15,12)$ lattices (Track B): our
  Theorem 2 predicts no PBB code there can exceed the CSS envelope at equal
  $n,k$; the search would test that prediction but has not been run.

**17:21 2026-08-18 sweep respawn** — the 3 nice-15 exp039 shards (35h, last cert Aug-17 13:38, ~1.7% CPU under box load) were SIGTERMed and respawned as 4 shards at nice 0: pids 12454/12456/12458/12460, `--ns 180 360 --stride 4`, logs results/partial_runs/exp039/shard_[0-3].log. Guard note: driver may print historical completion lines; certs land as parent_*.json mtimes — check mtimes, not log recency.
