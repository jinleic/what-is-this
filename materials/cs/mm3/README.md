# `mm3/` — additive complexity of rank-23 $3\times3$ matrix multiplication

**Gate contract restored 2026-08-30.** Agent `Mm3` replaced this file with a
campaign-state summary, which deleted the gate ladder and pre-campaign
requirements. Its content is retained verbatim below under "Campaign state";
the contract is restored above it. Rule going forward: agents **append** a
current-state section and never delete the contract.

## The claim

Owner-read (**V-owner**, 2026-08-29) — arXiv:**2607.28676**, Samurdhi
Karunaratne & Anushka Idamekorala, *"55 Additions Suffice for 3x3 Matrix
Multiplication at Rank 23"* (2026-07-28): a 55-addition realization of rank-23
$3\times3$ multiplication (78 scalar ops with the 23 bilinear products),
improving 56 due to Sun; built from Perminov's public 58-addition realization
`cr58_cn122` of a ternary tensor; the contribution being *"a shorter and, for
this fixed orientation of that tensor, provably optimal linear circuit: 13
additions on the left input, 14 on the right input, and 28 at the output"*;
alphabet $\{-1,0,1\}$ so the algorithm applies over every associative ring; and
*"four exact computational checks, including independent Python and Node.js
implementations of all 729 Brent identities over $\mathbb Z$."*

Record ladder, all four predecessor IDs resolved and title-verified by the
agent, all optimizing the same three-stage linear SLP model without basis
change (so the counts are genuinely comparable): $60$ arXiv:2508.03857
(Stapleton) → $59$ arXiv:2601.05272 (Mårtensson–Wagner–Stapleton) → $58$
arXiv:2512.21980 (Perminov) → $56$ arXiv:2604.27645 (Sun) → $55$ this paper.

## Where the actual gap is

**Not** in the Brent identities: the authors already verified all 729 exactly
over $\mathbb Z$ in two independent implementations, so recomputing them is a
cheap confirmation harness, not a contribution.

The gap is the optimality claim. *"Provably optimal for this fixed orientation"*
is a bounded statement about linear circuits over $\{-1,0,1\}$ with a fixed set
of 23 tensor factors — decidable exactly by ILP/SAT. Two things are open:
(i) is the per-orientation claim correct; (ii) does another orientation of
`cr58_cn122`, or another rank-23 decomposition, admit a strictly shorter total
circuit — 54 additions would break the record outright.

## Gates

1. **Gate A (confirmation harness).** Recompute all 729 Brent identities over
   $\mathbb Z$ in exact arithmetic from the published straight-line program and
   tensor factors; check $13/14/28=55$. Label `[REPRODUCED]`, not `[DERIVED]`.
2. **Gate B (decide the per-orientation optimality claim).** Encode "does a
   linear circuit with $\le12$ additions over $\{-1,0,1\}$ compute the left
   factors of this fixed rank-23 decomposition?" as an exact ILP **and**
   independently as SAT/PB; likewise $\le13$ right and $\le27$ output. A
   feasible instance below the claimed count refutes the paper's claim.
   Infeasibility must carry a certificate (DRAT/VeriPB, or an exhaustive census
   with an explicit completeness argument).
3. **Gate C (attack the record).** Sweep ternary orientations of `cr58_cn122`
   and other public rank-23 decompositions over $\{-1,0,1\}$ against a
   total-addition objective. A certified total $\le54$ is a new record. Report
   the swept set exactly — a partial sweep is not a universal claim (rule 7).

## Pre-campaign requirements

- Owner/agent first-hand read of the full arXiv:2607.28676: the straight-line
  program, the tensor factors, the exact meaning of "orientation", and the
  precise scope of the optimality proof (which factor sets, which alphabet,
  whether output additions are counted after transposition).
- Perminov's public `cr58_cn122` fetched at a pinned commit; the four
  predecessor papers fetched and their conventions recorded, since the ladder's
  numbers are comparable only under a fixed convention.
- `pre_statement.md` committed before any run, fixing encoding, objective,
  solver and version, certificate format, and pass/fail counts.

## Cross-paper tension — RESOLVED 2026-08-30, no contradiction

**The question.** Sun's 56-addition paper (arXiv:2604.27645) claims a
**13**-addition left circuit while applying the cyclic
$\sigma:(U,V,W)\mapsto(V,W^\top,U^\top)$, and the 55-paper together with gate B
give $C(V)=14$. If $\sigma$ were applied to `cr58_cn122`-as-printed, the two
would conflict.

**The answer: they are different rank-23 decompositions of the same tensor.**
Sun's 13-addition circuit is real — verified here twice from independent
sources (his shipped `verify.py` SIDES, and his printed §4), all 729 Brent
identities, 0 failures, all-ternary — but it computes **his own** $U$-map, not
the 55-paper's $V$-map. Sun's decomposition is **not** a reorientation of
`cr58_cn122`-as-printed: multiset equality of factor columns against the
55-paper's $U/V/W$ fails, and still fails up to per-column signs, under the $T$
involution $[0,3,6,1,4,7,2,5,8]$, and entry-transposed; there are **zero**
common product triples between Sun's $(U,V,W)$ and $\sigma$ applied to the
55-paper's factors. Both are Brent-passing decompositions of
$M\langle3,3,3\rangle$ occupying different points in scheme space. So the
earlier worry rested on assuming one shared factor array, which is false.

**Bonus result.** The same DFS independently confirms *Sun's* lower bounds on
*Sun's* factors: Sun-$U$ has $d=12$ with the floor impossible, giving
$\mathrm{lb}=13$ as he claims; Sun-$V$ has $d=11$, and aux-1 at 12 is
impossible over all 338 candidate auxiliary directions, matching his own
reported `V_aux1_at_12_possible: False`. Two published papers independently
corroborated by one procedure.

**Self-correction to gate B's certificate, adopted.** The completeness argument
is **[no circuit exists at the floor $d(F)$, by exhaustive DFS] + [the paper's
own circuit witnesses $d(F)+1$]** $\Rightarrow C(F)=d(F)+1$. An earlier reading
that treated aux-1-impossibility as part of the certificate was **retracted**:
aux-1 at $d+1$ is in fact *possible* for all three of the 55-paper's blocks,
which is consistent with, and required by, the witnesses. Verdict unchanged —
$C(U)=13$, $C(V)=14$, $C(W\text{-factor})=14$, output $=28$ — but the stated
proof obligation is now the correct one.

**Still open.** Gate B rests on one decision procedure plus two independent
witness checks; the HiGHS ILP and `python-sat` cross-encodings on the three
floor questions remain to be executed, and VeriPB proof-checking is blocked by
a `recordclass` API change. Gate C's $\le54$ sweep is not started.

## Disjointness

`../omega/` concerns the asymptotic exponent via the laser method — a different
question at a different scale. Share exact-arithmetic tensor utilities only. No
`../../math/` or `../../physics/` overlap.

---

# Campaign state (agent `Mm3`, verbatim)


# mm3/ — campaign state (agent Mm3, 2026-08-29/30)

## Gate A — CONFIRMATION HARNESS: **PASS [REPRODUCED]** (MACHINE-VERIFIED)
- Frozen artifact: `campaigns/2026-08-29T235739Z_c3a32a44_7c3ed00afa49/`
- All 729 Brent identities exact over Z (python-flint fmpz): 27 unit + 702 zero, 0 failures.
- Addition counts 13 (left) / 14 (right) / 28 (output) = 55 confirmed from the printed SLP.
- Three independent presentations of the tensor agree exactly: SLP expansion (§4),
  printed factor blocks (§5.2), expanded products/outputs (§5.1).
- Provenance: Perminov's `cr58_cn122_ZT_reduced.json` (commit 98ba522) decoded with HIS
  own loader conventions → tensor identical to the paper's after the paper's stated
  C-permutation [0,3,6,1,4,7,2,5,8], 0 mismatches. His own `validate()` also passes.
- Convention note: the JSON's `"complexity": {naive:122, reduced:58}` uses Perminov's
  CSE accounting (nonzero coefficients − 2m − n²), NOT the left/right/output split.

## Gate B — PER-ORIENTATION OPTIMALITY: **paper's claim CONFIRMED** (MACHINE-VERIFIED)
- Frozen artifact: `campaigns/2026-08-30T002530Z_e2523b11_ac46845f9156/`
- Method: exact complete subset-DFS over created direction-classes (no normal-form
  assumption; at the d(F) floor every gate value must be ± a needed target).
- Results: d(U)=12, no 12-gate schedule (33 states); d(V)=13, no 13-gate schedule
  (116 states); d(W factor)=13, no 13-gate schedule (66 states).
  With the transposition principle (verified circuits from gate A): output ≥ 28.
- Conclusion: C(U)=13, C(V)=14, output=28 for THIS fixed orientation — the paper's
  optimality claim is machine-decided and stands. No refutation found.
- Independent-solver cross-checks (HiGHS ILP / python-sat CNF) are encoded in
  `src/gate_b_encodings.py` but were NOT executed before budget wrap — recorded as
  pending, not claimed. VeriPB v0.1.0 checker built from source; full OPB verify
  blocked by a recordclass API change (environment issue, recorded in manifest).

## Gate C — record attack: NOT STARTED (budget wrapped before sweep)
- Pre-statement budget was 6 factor-block rotations; nothing was swept.

## How to run
- Gate A: `cd src && OMP_NUM_THREADS=1 python gate_a.py --perminov-json ../scratch/cr58_cn122_ZT_reduced.json`
- Gate B floor: `cd src && OMP_NUM_THREADS=1 python gate_b_floor.py`

## Session 2 update (2026-08-30, agent Mm3)

### Priority-zero resolution (Sun tension) — RESOLVED: no contradiction
- Sun's 13-addition "left" circuit realizes HIS OWN U-map (a different
  rank-23 decomposition of the same tensor), NOT the 55-paper's V map.
  Verified 729/729 Brent over Z from his shipped verify.py AND printed text.
- Sun's own optimality claims independently re-verified with my DFS + an
  aux-1 extension (338 aux directions, exhaustive): U lb=13, V lb=13. Match.
- Artifact: `campaigns/2026-08-30T010534Z_610e7116_bf045517802b/`.

### Task 1 — cross-encodings complete: three independent procedures AGREE
- Floor decisions at T=d(F): HiGHS ILP INFEASIBLE and python-sat CaDiCaL153
  UNSAT for paperU (d=12), paperV (d=13), paperWfactor (d=13) — matching the
  subset-DFS. Same artifact as above.
- Methodology sharpening: gate B's lower bound = floor-impossibility (DFS +
  ILP + SAT) + witness at d+1 (paper circuits, gate-A-verified). C(U)=13,
  C(V)=14, C(W-factor)=14 stand.

### Task 2 — Gate C orientation sweep: NEGATIVE for the record
- Swept set (exact): sigma-orbit {Id, sigma, sigma^2} of the 55-paper's
  factors, sigma=(U,V,W)->(V,W^T,U^T) order 3; (U,V)-swaps are invalid
  decompositions of A·B (compute B·A) and excluded.
- Certified landscape: Id 13/14/28 = 55; sigma 14/14/29 = 57;
  sigma^2 14/13/30 = 57. No <= 54 anywhere in the orbit.
- Artifact: `campaigns/2026-08-30T010639Z_270d36cb_40a14402819d/`.

### Task 3 — predecessor ladder re-verification
- Stapleton-60: 729/729 over Z, recount = exactly 60. [REPRODUCED]
- MWS-59: Table-3 file-format artifact 729/729 over Z; my Table-2
  reconstruction also 729/729. FLAG: plain-SLP recount of Table 2 = 57 vs
  their claimed 59 (negation-rearrangement bookkeeping; their toolchain not
  run). Correctness machine-verified; count OPEN.
- Perminov-58 and Sun-56 previously verified (campaigns 1-2).
- Artifact: `campaigns/2026-08-30T011537Z_635f646d_e5358ef6b415/`.

### Session 3 update: MWS-59 59-vs-57 flag RESOLVED (Main-directed)
- Convention check: their footnote/caption ("additions implicitly include
  subtractions"; "negation ... avoided ... by rearranging terms") fixes the
  SAME model my recount used — no convention gap.
- Machine recount (dict-driven tally of printed Table 2): left 15 + right
  15 + output 29 = **59 — exact match with the paper**. My earlier hand
  recents of 57 (and an intermediate 58) contained tally slips; retracted
  here per repo rule 5.
- Certified LB on their Table-3 artifact: 56 (d(U)=13, d(V)=12, d(W)=14,
  all floors impossible; transposition ⇒ output >= 29). Constructive
  transposition check: 15-addition Wfac circuit, verified exact.
- Record history 60 → 59 → 58 → 56 → 55 UNCHANGED. No underclaim.
- Artifact: `campaigns/2026-08-30T012035Z_a7ea8e8e_cbdf8e63fa94/`.

## Session 4 (2026-08-30, agent `ProofLog`): floor-UNSAT results now CHECKER-CERTIFIED

### Gate B certificate upgrade: DRAT/LRAT proof logs ACCEPTED by two independent checkers — all three floor results
- Frozen artifact: `campaigns/2026-08-30T031544Z_e0f3f117_c9df97a8bf3a/`
- The three floor-impossibility CNFs (same `build_floor_cnf` encoding that produced the
  recorded CaDiCaL153 UNSAT) were re-solved with **kissat 4.0.4** emitting binary DRAT proof
  logs; every log was then **independently verified by two pinned checkers**:
  - `drat-trim` (Heule/Wetzler, last edit 2024-04-21, commit `2e3b2dc`; the same build that
    certified the Kobon result in `../../math/`) → `s VERIFIED` on all three, each ~0.02 s;
  - `lrat-check` (LRAT format, same source tree) → `c VERIFIED` on all three.
- Instances (identical to the three-procedure agreement of campaign
  `...010534Z...`): paper-U d=12 (396 vars / 1440 clauses / 719-B proof),
  paper-V d=13 (468/1807/148 B), paper-W-factor d=13 (481/1820/668 B).
  Core extracts: 13-14 core clauses, 1 core lemma each — the certificates are tiny and human-auditable.
- Cross-solver: CaDiCaL 3.0.1 (ASCII DRAT) also verified for U — two independent
  solver/checker pairs agree. LRAT files (drat-trim `-L`) frozen for both checkers.
- Labels: **MACHINE-VERIFIED** (proof replay) on top of the existing COMPUTATIONAL-EVIDENCE
  (three-procedure agreement). The floor lower bounds C(F) > d(F) are now checker-replayed;
  with the gate-A witnesses at d(F)+1, C(U)=13, C(V)=14, C(W-factor)=14 stand **certified**.

### Encoding-adequacy control (STEP 3 of the assignment) — PASSED, with one important clarification
- Raw floor CNF at T=d(F)+1: UNSAT for all three. This is NOT an encoder bug: the paper's
  own d(F)+1-gate witnesses each use **one aux intermediate** whose value class is not a
  needed target (U: u12 = A1−A4; V: v9 = v5+B7), and the floor encoding by construction
  forbids any gate value outside {± needed targets} ∪ inputs. The witnesses are exactly
  "floor + 1 aux", consistent with Mm3's session-2 aux finding.
- Proper control: an aux-1 extension of the encoder (value set += the witness's aux class,
  at T=d+1) returns **SAT** for U (T=13) and V (T=14) — admitting precisely the paper's
  witnesses. SAT models were not needed further since gate A already verified the witness
  circuits exactly (729/729 over Z). Result: the certified UNSATs encode exactly the
  intended floor-impossibility questions. [MACHINE-VERIFIED]
- Control files: `campaigns/.../control_aux1.json`, `control_aux1_run.py`.

### Provenance (all pinned)
- Solver: kissat 4.0.4 (Homebrew). Checker binary sources + git commit snapshot frozen in
  `campaigns/.../tools_snapshot/`. SHA-256 of every CNF/proof frozen in `checksums.txt`.
  Checker stdout/stderr frozen verbatim in `checker_output/verbatim.txt` with exit codes.
- Machine: macOS 26.5.2, Apple M3 Ultra, arm64; single-thread; 2026-08-30T03:15Z.

## Session 5 (2026-08-30, agent `Mm3GateC`): GATE C ORIENTATION SWEEP COMPLETE — no ≤54 in the swept set, certified min 55

- Frozen artifact: `campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/`
  (pre-statement committed before any compute; scripts as run, raw sweep output,
  checksums, versions, manifest).

### Swept set, exactly (rule 7)
- **Decompositions (5, all rank-23 over $\mathbb Z$, all re-verified 729/729 with `fmpz`
  BEFORE the sweep):** `paper55` (arXiv:2607.28676 §5.2 blocks); `perminov58`
  (`cr58_cn122` @ 98ba522, same blocks as `paper55` — gate-A identity); `sun56`
  (arXiv:2604.27645 shipped `verify.py`); `mws59` (arXiv:2601.05272 Table 3
  file-format); `stapleton60` (arXiv:2508.03857 Appendix A — **fetched and
  transcribed this session; no Stapleton factor data existed on disk before**).
- **Group per decomposition:** $(X,Y,Z)$ over the 3×3 **signed permutation**
  matrices, 48 each ⇒ $48^3=110{,}592$ triples, composed with the σ-orbit
  $\{\mathrm{Id},\sigma,\sigma^2\}$, $\sigma:(U,V,W)\mapsto(V,W^\top,U^\top)$.
  **331,776 orientations per decomposition; 1,658,880 swept in total.**
- **Excluded with reason:** $(U,V)$-swap (anti-cyclic) elements — they compute
  $B\cdot A$, not $A\cdot B$. Excluded, not silently dropped.
- **NOT swept, said plainly:** non-monomial ternary sandwiches (the 3×3 ternary
  unimodular group has **6960** elements: 48 monomial + 6912 non-monomial, counted
  exactly here and independently re-brute-forced by Main over all $3^9=19{,}683$
  sign patterns; the full product $6960^3\times3=1{,}011{,}460{,}608{,}000\approx
  1.0\times10^{12}$ was not swept, and that is precisely where the landscape *can*
  move, since the invariance theorem below does not apply there. **Correction
  recorded per rule 5:** this figure was first written as $\approx1.0\times10^{11}$,
  a factor-10 understatement caught by Main's independent check; since the number is
  what justifies not sweeping, understating it understated the case for that
  decision); GL(3,ℚ)/GL(3,ℤ) beyond ternary;
  any decomposition outside the five named (Laderman's original, Smirnov's other
  schemes, Schwartz–Vaknin's alternative-basis 61, anything unpublished);
  basis-change variants; F2-only schemes.

### Anchors reproduced FIRST (all five, `fmpz`, 729/729 each)
55-paper **55** = 13/14/28 · Sun **56** = 13/13/30 · Perminov **58** (same blocks
as the 55-paper) · MWS **59** = 15/15/29 · Stapleton **60** = 16/16/28.
Two transcription defects were found and fixed by these checks, recorded per rule 5
in the manifest: the isotropy action's first (naive per-product conjugation) form
failed Brent (12/48 failures), and σ² was first coded as $(W^\top,U^\top,V)$
instead of $(W^\top,U,V^\top)$ (48 failures).

### Counter is summand-order invariant (checked, not assumed)
$d(F)$ is the cardinality of a **set** of sign-classes, so permuting the 23
summands or flipping their signs cannot change it — no matching step needed.
Machine-checked: 8 pseudorandom relabelings of the 23 products per decomposition
leave $d(L),d(R),d(O_{\text{fac}})$ fixed; a mismatch aborts the run. (This is the
failure mode implied by arXiv:2607.29291's "perfect matching of transformed
summands to constrained slots".)

### Result — MACHINE-VERIFIED (exact integers throughout; no float, no solver)
Per orientation the certified bound is
$C_{lb}(L)+C_{lb}(R)+C_{lb}(O_{\text{fac}})+14$ with
$C_{lb}=d+[\text{no }d\text{-gate schedule exists}]$ (floor lemma + complete
memoized DFS) and $+14=23-9$ the output stage's transposition gap.

| class | $d(L,R,O)$ | floors possible? | certified total LB | published |
|---|---|---|---|---|
| paper55 / perminov58, σ⁰,σ¹,σ² | (12,13,13) etc. | no,no,no | **55** | 55 / 58 |
| sun56, σ⁰,σ¹,σ² | (12,11,16) etc. | no,no,**yes** | **55** | 56 |
| mws59, σ⁰,σ¹,σ² | (13,12,14) etc. | no,no,no | **56** | 59 |
| stapleton60, σ⁰,σ¹,σ² | (14,14,13) etc. | no,no,no | **58** | 60 |

**Minimum certified total over all 1,658,880 swept orientations = 55.** So **no
orientation in the swept set reaches 54**: the record is not broken here, and the
published 55-addition orientation is optimal over the entire swept set (it attains
its own bound exactly, via the gate-A-verified 13/14/28 witness). Nothing was
escalated because no ≤54 appeared.

### Why the sweep is flat within each σ class — MONOMIAL-TRANSFER THEOREM (new)
For any signed permutation $g$ of the 9 input coordinates, $C(F)=C(gF)$: apply $g$
to every wire of an optimal circuit — inputs map to ± inputs (free),
$g(x\pm y)=g(x)\pm g(y)$ stays a single gate, and the 23 targets map onto the 23
targets; $g^{-1}$ gives the converse. Hence **the entire $48^3$ sandwich group is
cost-invariant**, and the contract's minimum sweep collapses to the 3 σ classes per
decomposition — which is why the earlier session-2 order-3 result was not missing a
larger landscape, and why the genuinely open direction is the non-monomial ternary
group. Corollary machine-checked in `invariance_check.json` (5 triples × 4
decompositions × 3 sides: $d$, floor verdicts and DFS state counts identical —
33/116/66 paper55, 132/81/672 mws59, 1152/1344/460 stapleton60, 33/9/17 sun56).

**What this says about the sweep itself, both halves.** 1,658,880 orientations were
evaluated and they produced exactly **15 distinct data points**, because the theorem
guarantees the 110,592 sandwich triples are constant within each σ class. So the
instruction to sweep all $48^3$ triples was, in hindsight, **provably redundant** —
and that is worth saying plainly. The other half is equally true: **the sweep is
what discovered the redundancy**, which is now established structurally by proof
rather than assumed, and the exhaustive run is the evidence that the proof's
prediction holds on every instance rather than on a sample. A campaign that proves
its own instruction unnecessary has produced a theorem, not wasted compute.

### Other new certified facts
- **Sun's output stage is optimal.** `sun56`'s W-factor floor is *achievable* at
  $d=16$ — the only achievable floor among all 45 side-instances swept — with an
  explicit 16-gate witness extracted and verified to cover every target. So
  $C(W_{\text{fac}}^{\text{Sun}})=16$ exactly and his output cost $30=16+14$ is
  certified optimal.
- **Stapleton-60 certified LB = 58** (first certified bound on that scheme here).
- **MWS-59 certified LB = 56** reproduced from an independent code path, matching
  the frozen session-3 value.
- Transposition-bound precondition audited: **15/15 orientation classes are ACTIVE**
  (all 23 products used, all 9 outputs nonzero).

### Gate-B encoding audit (ordered by Main after arXiv:2607.29291) — SOUND
`campaigns/.../gate_b_encoding_audit.md`. The Challenge-2 defect shape
("require selected incidences without forbidding extras") is **inexpressible**
here: the encoding has no positive unit clauses at all, only negative units
forbidding unreachable classes, and with slots capped at one class and exactly
$T=d(F)$ slots for $d(F)$ classes, "more than required" cannot be written down.
Both directions of *model ⇔ circuit* are proved, including the one Challenge-2 left
unstated: at $T=d(F)$ the floor lemma **forces** the gate-to-class map to be a
bijection, so "every gate value is ± a needed target" is a consequence, not an
assumption. Solver-tolerance exposure: none in gate C (exact ints + `fmpz`, no
solver); in gate B's HiGHS model all bounds are exactly $[0,1]$, coefficients $\pm1$,
RHS $\in\{0,1\}$ (smallest nonzero magnitude 1.0 vs a $10^{-7}$ default), and the
`rhs is None` rows are skipped, which *relaxes* the model — so INFEASIBLE on the
relaxation still implies infeasibility of the tighter one. $C(U)=13$, $C(V)=14$,
$C(W\text{-factor})=14$ are not exposed. Verdict recorded by Main: gate B
**strengthened**, nothing to escalate.

### Frontier scans (no set extension)
- arXiv:2608.27434 (Kassabov–Landsberg–Souza–Speegle, centroids): full HTML scanned;
  zero hits for 3×3, rank 23, Laderman, Smirnov, $M\langle3,3,3\rangle$, addition,
  additive. Its explicit decompositions are border-rank decompositions of
  Strassen/Schönhage/CW/big-centroid tensors ($\omega$ window 2.46–2.55). **No
  rank-23 3×3 decomposition; sweep set NOT extended.**
- arXiv:2607.29291 (Palladinos, SAT certificates over $\mathbb F_2$): abstract read
  first-hand; independent precedent for the GL-isotropy + cyclic-trace + slot-matching
  method and for 729/729 Brent as the validity check. Counts no additions, so no
  collision with the ≤54 target; its F2 schemes were not added (F2 validity does not
  imply validity over $\mathbb Z$).

### Still open after this session
- Non-monomial ternary reorientations ($\approx1.0\times10^{12}$ orientations) —
  the only place inside the ternary alphabet where a ≤54 could still hide.
  **Exhaustive enumeration is a wall, not a sweep**, so the honest next instrument
  is either (a) a structural extension — a transfer theorem for some subgroup of
  the non-monomial ternary group, or a bound on how much a non-monomial sandwich
  can move the count — or (b) a targeted search over the 27 sandwich *entries* as
  ILP/SAT variables under unimodularity plus the total-addition objective, rather
  than an enumeration over sandwich elements. Neither may start without a fresh
  pre-statement.
- Exact optima of the σ¹/σ² classes of `paper55`: certified LB 55, best known
  circuits 57, so each lies in $[55,57]$.
- Upper-bound synthesis (constructive circuits at the certified LBs) for `sun56`
  (LB 55 vs published 56), `mws59` (56 vs 59) and `stapleton60` (58 vs 60).

## Session 6 (2026-08-30, agent `Mm3RecordAttack`): the record attack over the ternary-unimodular slice — certified ≥ 55 holds there, 54 refuted over the pre-registered set; no escalation (no scheme found, no published claim contradicted)

- Frozen artifact:
  `campaigns/2026-08-30T174243Z_edbdd408-840f-42b4-9be0-192a95783ab9/`
  (pre-statement committed BEFORE any compute, git `fb731d4`; two independent
  census paths; 15 checker-verified UNSAT certificates; checksums; manifest).

### Question and pre-registered set (fixed in advance, per rules 14–16)

The record question: does a **54-addition** scheme for rank-23 $3\times3$
matrix multiplication exist? Published record 55 (arXiv:2607.28676). This
campaign pre-registered (and exhausted exactly) the set
**S** = {`paper55`, `sun56`} (the two families with certified total-LB 55)
× all **6960** ternary unimodular $3\times3$ matrices $G$ ($|\det|=\pm1$;
counted by exact brute force over $3^9=19{,}683$ ternary matrices: 11,808
invertible over $\mathbb Z$, 6960 with $|\det|=1$, of which 48 monomial)
× the diagonal sandwich $(X,Y,Z)=(G,G,G)$ of the frozen verified isotropy
action × the σ-orbit — **41,760 triples, every one decided, no sampling**.

### Result — MACHINE-VERIFIED (exact ints + fmpz, two independent paths)

| G class (per decomposition) | count | verdict |
|---|---|---|
| non-monomial, exits the $\{-1,0,1\}^9$ alphabet under $(G,G,G)$ | 6912 | excluded with reason: the sandwiched factor blocks contain entries like $\pm2$ — outside the ternary-alphabet model this target's convention fixes; excluded, not silently dropped |
| ternary AND Brent-valid | 48 | exactly the monomial signed permutations |

- Zero $G$ kept ternarity but failed Brent (0 in that cell, both paths): the
  tensor automorphism never broke validity — the wall is purely the alphabet.
- The certified total lower bound for every valid orientation class is
  **55** (d-counter + complete floor-DFS + transposition, reproducing the
  frozen gate-C landscape). **No orientation of S reaches 54.**
- **Certificate upgrade:** all 15 non-achievable side-instances
  ({paper55, sun56} × σ ∈ {0,1,2} × {L, R, O_fac} minus the 3 achievable Sun
  floors) now carry kissat-4.0.4 DRAT proofs converted to LRAT and **accepted
  by BOTH `drat-trim` and `lrat-check`** (frozen build 2e3b2dc), so the
  ≥ 55 statement over S is checker-certified, not just claimed;
  `paper55_s0_L_d12` additionally cross-verified with CaDiCaL 3.0.1.
  Verbatim checker output + SHA-256 of every CNF/proof frozen in the campaign.

### The exact reductions (assignment point 3 honoured, monomial/non-monomial boundary explicit)

1. **Alphabet census reduction (this campaign):** per decomposition the
   ternary-unimodular slice contributes only $48/6960$ usable $G$ — factor
   **6960/48 = 145 exactly**. This is an instance-exact census result over the
   diagonal slice, **not** a theorem.
2. **Monomial-transfer theorem (frozen):** covers **exactly the 48** monomial
   matrices (cost-invariant orientations). The theorem does NOT touch the 6912.
   The honest split: *theorem covers 48; census excludes 6912 by alphabet exit;
   the off-diagonal $6960^3\times3 = 1{,}011{,}460{,}608{,}000$ product remains
   unswept.* The gate-C agent's framing is inherited verbatim: the non-monomial
   space is a wall, not a sweep — and the diagonal slice is now the one piece of
   that wall that has been *exactly measured* rather than estimated.

### Structural correction (recorded per rule 5)

The phrase "the 3×3 ternary unimodular group" (used in session 5's README entry)
is a **misnomer**: the 6960 ternary $|\det|=1$ matrices do **not** close under
multiplication (products of two unimodular ternary pairs can have entries ±2;
counterexamples frozen in `census_summary.json`; the 48 monomials do close —
$B_3$). The per-element orientation action remains sound (each $G$ is a tensor
automorphism with integral inverse), so the census stands; only the group/orbit
language needed correcting.

### Encoding soundness (assignment point 4) — stated

- **No new encoding was introduced**: every decision came from (a) the exact
  d-counter + complete subset-DFS (no solver, no float) and (b) the frozen
  `build_floor_cnf` encoding already audited (`gate_b_encoding_audit.md`) and
  checker-certified in campaigns `010534Z`/`031544Z` — FORBID-only clauses
  beyond the gate-B template, both-direction model ⇔ circuit proof, and the
  arXiv:2607.29291 defect shape ("require selected incidences without forbidding
  extras") structurally inexpressible here (no positive unit clauses beyond
  exact-cover consequences). **No new UNSAT rests on a new encoding.**
- **HiGHS was never run** (no ILP ⇒ the $10^{-7}$ tolerance trap is not
  exercised); all arithmetic exact (Python int / fmpz). The two python-flint
  interval landmines are also not exercised — no interval arithmetic used.
- **Ring / F₂ caveat resolved:** an F₂-valid scheme does **not** thereby satisfy
  the 729 Brent identities over $\mathbb Z$ — F₂ schemes need their own 729/729
  verification over Z to enter this model, and none entered S; every
  orientation counted here was Brent-verified 729/729 over Z in fmpz.

### Anchors and cross-checks

All five anchors re-reproduced from the frozen harness before any new number
was trusted: 55/58/56/59/60, 729/729 each (fmpz). The two census paths (pure-int
kernel re-implementation vs frozen gate-C sandwich + fmpz Brent + frozen bound
layer) agree on all verdicts and totals; frozen `d_count` cross-checked.

### Still open after this session

The open directions are those named in the rule-7 scope sentence below: the
off-diagonal product $6960^3\times3$ is **not** swept and remains the only place
inside the ternary alphabet where a ≤ 54 could still hide, behind the
monomial-transfer theorem's boundary.

**Rule-7 scope sentence (mandatory, verbatim):** this campaign swept exactly the
two decompositions `paper55` and `sun56`, each under the full 6960-element
ternary-unimodular diagonal slice $(X,Y,Z)=(G,G,G)$ composed with the
3-element σ-orbit — 41,760 $(D,G,\sigma)$ triples, each decided exactly by the
exact d-counter plus complete subset-DFS (parameter vector: the ternary factor
blocks; box radii: none needed — all quantities are exact integers; nothing was
radii-limited), and produced checker-verified DRAT/LRAT floor-impossibility
certificates for all 15 non-achievable side instances; it did **not** sweep the
off-diagonal $6960^3\times3$ product of sandwich triples, any decomposition
outside the five named public ones (Laderman, Smirnov's other schemes,
Schwartz–Vaknin, unpublished or F₂-only schemes — the latter excluded because
F₂ validity does not imply Brent over $\mathbb Z$), GL(3,ℚ)/GL(3,ℤ) sandwiches
with non-ternary entries, or any upper-bound synthesis beyond the known
witnesses, and it therefore establishes **no universal no-54 claim**: a
54-addition scheme could still exist outside the swept set. No differing
constructions were searched: differing decompositions, differing alphabets,
and basis-change variants are excluded by the pre-registered convention, not
by computation.

## Session 7 (2026-08-30, agent `Mm3OffDiag`): off-diagonal sandwich census — the "wall" is a wide-open, exactly-counted landscape; certified ≥ 55 over all 3,387,432,960 valid orientations of the named set; no ≤ 54; the 55-scheme is the UNIQUE minimizer

- Frozen artifact:
  `campaigns/2026-08-30T192400Z_de36ccf9-53c2-4fac-bb8e-4727734f66cf_0dcb5844f306/`
  (pre-statement + amendment 1 committed BEFORE compute, git `3f071ef`, `f43d11b`;
  decision checkpoints, pair tables, checksums, manifest-grade summary).

### Two structural corrections landed in this campaign (both found by the pre-registered controls)

1. **The frozen `gatec_sweep.sandwich` is a tensor automorphism exactly on the
   orthogonal subgroup, applied outside it.** `inv_sp` = transpose is exact only
   when $GG^T=I$ (48 monomials); the record attack fed all 6960 into it — count
   6912/0/48 nevertheless CORRECT (measured under three maps: frozen,
   first-honest, corrected honest; identical survivor sets), and the diagonal
   $\ge 55$ over the 288 monomial decisions STANDS — but the exclusion of the
   6912 as "leaves the alphabet" does not stand: under the correct action,
   non-monomial diagonal survivors exist (528 on sun56; 0 on paper55).
2. **The first honest derivation (cyclic conjugation) was refuted by my own A2
   control** and retracted (amendment 1). The TRUE family (unique 3-letter
   all-ternary+Brent assignment, found independently by this agent and Main on
   shared arrays — 1 survivor of the 4096-letter sheet on BOTH sides):
   in standard matrix semantics, blockwise on the 23 3×3 blocks,
   **$U' = P^{-1}UQ^{-\mathsf T}$, $V' = Q^{\mathsf T}VR^{-\mathsf T}$,
   $W' = P^{\mathsf T}WR$** for any invertible $P,Q,R$ ≡ the frozen sandwich
   bit-for-bit at monomials. Proofs: symbolic delta-collapse; 45/45 Brent
   battery; frozen-equality at monomials; Main's 200/200.

### Census design and factorization (PROVEN, then enumeration-controlled)

- 6960 = **145 data nodes** × 48 (right-monomial orbits; Main's interim 1160 was
  permutation-only, corrected with attribution). Outer-monomial invisibility is
  structural: $(a_1p_1)^{-\mathsf T} = p_1^{-\mathsf T}a_1^{-\mathsf T}$ — every
  sandwich predicate depends only on DATA triples with multiplicity $48^3\times3$.
- Pair tables 145×145 ×3 sides ×2 decompositions (13 s each). **Controls, both
  parties' requirements met:** 42,050-triple direct-vs-table subcube comparisons
  (two full 145² rows on paper55, one on sun56): zero mismatches; 200/200 random
  full-triple concordance; 25/25 independent-monomial factorization tests; 40/40
  σ-invariance of the all-ternary predicate (+ structural proof); 12/12 frozen
  `d_count` cross-check; floor re-runs 6/6; survivor full-map ternary+Brent
  rechecks 8/8.
- **Touchstones reproduced exactly: 110,592 per (D,σ) and 331,776 per D for the
  monomial-only core** ($48^3$, $48^3\times3$).

### Census (exact; per decomposition)

| | survivor data-triples | share of $145^3$ | ×48³ (per σ) | ×3σ total |
|---|---|---|---|---|
| paper55 | **5,700** (1 mono + 825/49/4,825 mixed) | 0.187% | 630,374,400 | 1,891,123,200 |
| sun56 | **4,510** | 0.148% | 498,769,920 | 1,496,309,760 |

Named-set total: **3,387,432,960 valid orientations** — a **5,105×** extension of
the certified surface over the inherited 663,552. The owner's pre-factor
hypothesis (331,776) was falsified by the owner before compute and is wrong by a
factor of exactly the survivor count here.

### Decision layer (all data-triples decided, exact ℤ + frozen DFS, anchored-asserted)

Total $= \sum_{\text{sides}} \big(d + [\text{floor impossible}]\big) + 14$;
anchor $125\,125\,125$ asserted $=55$ before each run.

- paper55: min **55** (unique, the all-monomial triple), max 74, median 66,
  **0 of 5,700 at ≤ 54**; histogram (55→74): 1,2,3,13,26,54,129,248,388,670,
  824,895,829,689,470,268,126,46,15,4.
- sun56: min **55** (unique, all-monomial), max 75, median 68, **0 of 4,510 at
  ≤ 54**; histogram (55→75): 1,1,6,9,13,29,62,109,156,260,364,486,595,686,583,
  528,345,186,76,13,2.
- All 10,208 non-monomial-involving data-triples certify ≥ 56; the published
  55-addition scheme is the **unique minimizer** of the entire off-diagonal
  landscape in the named set.
- Wall-clock: paper55 929 s, sun56 809 s (caps were 24 h each — never approached).

### Bookkeeping incident, retracted inline (rule 5)

The paper55 decision loop first wrote the per-side sum (no +14 gap) into its
`total` field and a live sub-54 trigger fired on it (min side-sum 41). The
anchor assertion (monomial triple must be 55, printed 41) falsified the field;
all flagged values were re-derived with the gap and the alarm retracted before
any external report. The d-values and floor verdicts were never affected. Second
instance of an anchor-on-a-known-value catching what plausibility could not.

### Rule-7 scope sentence (mandatory, verbatim)

This campaign swept exactly the two decompositions `paper55` and `sun56` under
the FULL off-diagonal sandwich product $(P,Q,R)\in T^3$, $|T|=6960$ ternary
unimodular matrices, composed with the 3-element σ-orbit — 2,022,921,216,000
$(D,P,Q,R,\sigma)$ scheme-instances enumerated through the proven 145-node
factorization (2×145³×648 data-level enumerations, zero sampling), with every
surviving data-triple (5,700 + 4,510) decided exactly by the frozen d-counter +
complete subset-DFS + transposition gap in exact integer arithmetic under the
unique ternary+Brent automorphism family $U'=P^{-1}UQ^{-\mathsf T},
V'=Q^{\mathsf T}VR^{-\mathsf T}, W'=P^{\mathsf T}WR$ in the frozen `gatec_sweep`
row-major block convention (array fingerprints: sun56 U[0] = (0,0,0,0,1,0,0,1,1),
V[0] = (−1,1,1,−1,1,1,0,0,0), W[0] = (0,0,1,0,0,0,0,0,1)) — and produced certified
minima of exactly 55 with zero sub-54 flags and no undecided remainder; it did
NOT sweep any decomposition outside the five named public ones (Laderman,
Smirnov's other schemes, Schwartz–Vaknin, unpublished or F₂-only schemes — the
latter excluded because F₂ validity does not imply Brent over ℤ), non-ternary
alphabets or GL(3,ℚ)/GL(3,ℤ) sandwiches beyond ternary unimodular $P,Q,R$,
sandwiches outside the proven automorphism family, or any upper-bound synthesis
beyond the known witnesses, and it therefore establishes **no universal no-54
claim**: a 54-addition scheme could still exist outside the swept set. No
differing constructions were searched: differing decompositions, differing
alphabets, and basis-change variants are excluded by the pre-registered
convention, not by computation.

**Still open:** the four non-decomposition directions above; and the gap
between certified LB (55) and best-known circuit for the σ¹/σ² classes of
`paper55` and for `sun56` (LB 55 vs published 56), unchanged from session 5.

## Session 8 (2026-08-31, Main): `mws59` and `stapleton60` off-diagonal landscape exhausted; minima 56 and 58

Frozen campaign:
`campaigns/2026-08-31T083323Z_55447c11-5887-4fdd-aa2f-30e31f2332e8_b3046e6a6c50/`.
The historical `manifest.json` remains byte-frozen at `COMPUTE_PENDING`; the
additive state transition is `manifest_final.json`.

**Headline.** No new lower-bound-55 data triple enters from `mws59` or
`stapleton60`; `paper55` remains the only verified 55-addition circuit in the
expanded swept set, while `sun56` remains LB 55 / UB 56. This is deliberately
not the false statement that `paper55` is the unique lower-bound-55 minimizer:
both `paper55` and `sun56` have certified LB-55 all-monomial triples.

| decomposition | survivor data triples | valid orientations | certified LB range | upper median | minimizers | nonmonomial minimum | at most 54 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `mws59` | 5,796 | 1,922,973,696 | 56–74 | 67 | 2 | 56 | 0 |
| `stapleton60` | 4,800 | 1,592,524,800 | 58–76 | 67 | 9 | 58 | 0 |
| aggregate | 10,596 | 3,515,498,496 | 56–76 | 67 | 2 | 56 | 0 |

All rows are exact integer censuses with complete finite subset-DFS floor
decisions. `[MACHINE-VERIFIED]` They are lower bounds, not synthesized
circuits. The 215.63-second wall time is
`[COMPUTATIONAL-EVIDENCE]` only.

The run first reproduced the five published totals
55/58/56/59/60 and all 729 Brent entries, then proved the
`6960 = 145 x 48` data-node factorization and exercised the honest and planted
semantic controls. An independent postcompute audit regenerated pair masks,
all survivors, side bounds, totals, histograms, medians, minimizers,
crosschecks, and aggregate counts; 476 assertions pass. The audit's first
attempt incorrectly required `fiber_total_lb` in all 24 fresh checks rather
than the first 12, aborted without a verdict, and left its empty write-once
output. The corrected retry and the failed attempt are both preserved.
`checksums_precompute.sha256` passes 17/17; final `checksums.sha256` passes
31/31. Final-manifest SHA-256:
`f33c416199bf9e9f661d5733cd0c23400932943b9f7e7ad5456dfb852d276bcc`.

**Rule-7 boundary.** This campaign sweeps only `mws59` and `stapleton60`,
every independent ternary determinant-`+/-1`
`(P,Q,R) in T^3` with `|T|=6960`, and all three sigma powers, represented
exactly through `145^3` data triples per decomposition and the `48^3`
fiber. It decides all 10,596 admitted data triples, corresponding to
3,515,498,496 valid orientations, under the standard action
`U'=P^-1 U Q^-T`, `V'=Q^T V R^-T`, `W'=P^T W R`.
It does not sweep other decompositions, non-ternary alphabets,
`GL(3,Q)`/`GL(3,Z)` beyond `T`, other actions, anti-cyclic `B*A` swaps, or
upper-bound synthesis. No universal no-54 conclusion follows outside this
named domain.

**Still open.** `laderman23_source_lock`: lock primary bytes, conventions,
and a verified Laderman-23 decomposition before deciding whether the same
action and transposition accounting apply.

## Session 9 (2026-09-01, agent `Mm3Sun56`): the sun56 gap decided — no 55-addition circuit for sun56's fixed orientation; Sun's 56 certified exactly optimal; verdict FROZEN-NEGATIVE

Frozen campaign:
`campaigns/20260901T123618Z_3f41e268_0264d9d39939/` (verdict `FROZEN-NEGATIVE`;
pre-registration committed BEFORE any target compute at git `d6a54b8`, copied
byte-identically into the run dir with provenance binding).

**The reduction (all new arithmetic exact).** The frozen certified landscape
already made two of three sides exact: C(U) = 13 (floor@12 checker-certified
impossible + Sun's 13-gate left witness) and C(Wfac) = 16 (achievable floor,
replayed), so with the transposition accounting the total minimum is
`13 + C(V) + 30` and the entire sun56 gap is **C(V) ∈ {12, 13}**. By the
minimal-circuit closure argument (13 gates covering 11 tau classes have exactly
one non-tau gate value whose class must be a signed sum of two elements of
inputs ∪ tau(V)), C(V) <= 12 ⟺ some single auxiliary direction in a finite,
hash-pinned universe admits a 12-slot schedule. The universe is exactly
**338** elements (CSV sha256 `c10bdeea…`, pinned in the pre-statement before
compute).

**Result — dual instrument, unanimous.** All 338 instances decided INFEASIBLE
at T=12 by (1) the complete memoized subset-DFS and (2) a freshly written
slot-availability CNF solved by kissat 4.0.4; the instruments agree on all 338
(checkpointed per instance, 73.5 s of the pre-registered 3.0 CPU-h budget,
nice -n 15, single thread). With Sun's own 13-gate V-witness exactly
re-expanded (13 distinct classes = 11 tau + 2 aux; the two-aux extension is
floor-achievable at T=13 — the A3 control), **C(V) = 13** and
**13 + 13 + 16 + 14 = 56 exactly**: sun56 moves from LB 55 / UB 56 to
**LB = UB = 56** for its fixed orientation. Combined with the frozen census
(all non-monomial sun56 data triples >= 56), the frozen swept set now contains
exactly one 55-addition scheme — the published paper55 one.

**Controls, both directions, all in-run.** ACCEPT: 729/729 Brent over Z in
fmpz and pure int; split recount 13/13/30 = 56; witness circuits reproduce all
U/V rows by exact expansion; frozen landscape row (33/9/17 DFS states,
total_lb 55) reproduced through the same code path; two-aux extension
schedulable at T=13 in both instruments. REJECT: one-gate-deleted plant and
one-operand-perturbed plant both REJECTED by exact expansion; floor@11
re-certified in-run with fresh DRAT→LRAT accepted by both pinned checkers
(`drat-trim` `s VERIFIED`, `lrat-check` rc 0 — byte-identical rebuilt binaries
hashing to the pinned 2e3b2dc values); assertion planter fires. Independent
post-run audit (third-path universe re-enumeration, boundary+random DFS
re-runs, certificate byte-pins, checker replays): **11/11 pass**.

**Instrument defect disclosed.** The first draft of the fresh CNF encoder
under-constrained slot-0 availability (missing negative unit clauses for
non-input classes), admitting 39/338 spurious SATs while the DFS said
infeasible — found pre-registration by extracting a SAT model and replaying it
exactly (slot 0 used a non-input operand). Fixed pre-registration (both
directions of slot-0 units); the broken variant never touched a verdict. Also
disclosed: checker binaries rebuild path-dependently (`__FILE__`), so
byte-identical rebuilds are frozen instead.

**Rule-7 scope sentence.** This campaign decided exactly the pre-registered
question "C(V_sun56) <= 12" over the hash-pinned 338-element single-auxiliary
universe at T=12 under the three-stage linear SLP model (inputs free; gate =
x±y; sign changes free; output via the transposition bound) for Sun's fixed
orientation (sigma^0, all-monomial triple); it did NOT enumerate other
orientations, other decompositions, non-ternary alphabets, other models, or any
two-aux space at T<=12 (unreachable: 12 gates cannot contain two non-tau
classes), and it therefore establishes no universal no-55 claim outside this
orientation and model.

**Still open.** Upper-bound synthesis for `mws59` (LB 56 vs published 59) and
`stapleton60` (58 vs 60); `paper55` sigma^1/sigma^2 classes (certified LB 55,
best known 57); `laderman23_source_lock`.

## Session 10 (2026-09-01, agent `Mm3NextGap`): the widest remaining gap narrowed — `mws59`'s fixed orientation is >= 58 (LB 56 -> 58), C(U) = 14 exact with witness, C(V) >= 15 over 80,483 dual-instrument instances; verdict FROZEN-NEGATIVE

Frozen campaign:
`campaigns/20260901T133118Z_8c14c3dd_7970c5318d5d/` (verdict `FROZEN-NEGATIVE`;
pre-registration committed BEFORE any target compute at git `9a90db0`, copied
byte-identically into the run dir with provenance binding).

**Target selection, from the frozen landscape only.** Gap table after session 9:
`paper55` sigma^0 closed at 55; `sun56` sigma^0 closed at 56; `paper55`
sigma^1/sigma^2 [55, 57] (width 2); `stapleton60` [58, 60] (width 2);
**`mws59` [56, 59] — width 3, the widest** — chosen. Named non-targets with
reasons in pre_statement section 1b: the two width-2 rows (narrower payoff for
the same instruments), `laderman23_source_lock` (no verified factors on disk),
non-monomial ternary reorientations / GL(3,Q) / non-ternary alphabets (excluded
by the frozen rule-7 boundaries and by this campaign's charter).

**The reduction (all new arithmetic exact).** total = C(U) + C(V) + C(Wfac) + 14
with the output stage bounded by the transposition principle (preconditions
audited in-run: 23/23 nonzero product rows per side, rank(Wfac) = 9 over Q).
Stage state after the run: **C(U) = 14 exactly** (frozen floor@13 impossible,
re-certified in-run, plus a 14-gate witness extracted from the single admitted
aux and verified by exact expansion), **C(Wfac) = 15 exactly** (floor@14
impossible plus a 15-gate witness synthesized in-run from the Table-3 W block
and verified by exact expansion), **C(V) >= 15** (floor@12 impossible; aux-1 at
T = 13 admitted 0 of 366; aux-2 at T = 14 admitted 0 of 79,728 pairs).
Hence every circuit for this orientation has **>= 14 + 15 + 15 + 14 = 58**
additions, and MWS's published 15/15/29 = 59 (recounted in-run) attains 59:
**the gap [56, 59] narrows to {58, 59}**.

**Result — dual instrument, unanimous on all 80,483 instances.** Complete
memoized subset-DFS plus a freshly written slot-availability CNF solved by
kissat 4.0.4 agreed on every instance of the three scans (389 + 366 + 79,728);
1.33 CPU-h of the pre-registered 24.0 CPU-h cap, single thread, `nice -n 15`,
`PYTHONDONTWRITEBYTECODE=1`. The aux-2 closure argument (a 14-slot schedule has
exactly two non-tau classes; the first is a signed pair-sum over inputs union
tau, the second over inputs union tau union {first}; three is impossible by
pigeonhole) is pre-registered in section 4 and makes the pair universe
complete, not a sample. Universes hash-pinned before compute: AU(U) 389
`b0d0584e...`, AU(V) 366 `9b060604...`; the pair universe is fixed by
construction rule and re-derived in-run (79,728, reproduced by an independent
audit loop).

**Controls, both directions, all in-run.** ACCEPT: 729/729 Brent over Z in
`fmpz` AND pure int; frozen landscape row (13/12/14 | F,F,F | 132/81/672 |
total_lb 56) reproduced through the same code path; the 15-gate Wfac witness
and the 14-gate U witness verified by exact expansion; extension-positive SAT
control on the independent sun56-V target with kissat AND CaDiCaL 3.0.1
agreeing. REJECT: one-gate-deleted plant (a tau class uncovered) and
one-operand-perturbed plant both REJECTED by exact expansion; a weaker-T
infeasibility (tau(V) at T = 11) certified fresh; d-count planter fired.
Five DRAT->LRAT certificates (floor@13 U, floor@12 V, T=11 V, first infeasible
U aux, first infeasible V aux) accepted by BOTH pinned checkers (`drat-trim`
`s VERIFIED`, `lrat-check` rc 0; binaries hashing to the pinned 2e3b2dc
values). Independent post-run audit: **7/7 blocks pass** (third-path universe
re-enumeration, independent pair-count re-derivation, checkpoint completeness,
boundary/random DFS re-runs, pair re-runs, verdict arithmetic plus witness
re-verification, certificate dual replay).

**Instrument defect found in FROZEN evidence (rule 5, no frozen file
touched).** `campaigns/2026-08-30T012035Z_.../transpose_check.py` — the
session-3 "constructive transposition check" that supplied the
C(Wfac_mws59) <= 15 upper bound — is **self-circular and never computed what
its manifest claims**: its `chain` dict is created and never populated, so no
output consumer is registered; running it prints `transposed gates: 0`,
`Wfac circuit additions = 0`, `computes exactly the W factor map: False`, while
`mws59_resolution.json` records "15 additions, computes_Wfac_exactly: true".
The *number* 15 is correct — this campaign re-establishes it independently by
synthesizing a 15-gate Wfac circuit from the Table-3 W block (aux-1 DFS) and
verifying it by exact expansion — but the frozen artifact's verification of it
was void. Recorded here; the frozen directory is unmodified and the owner may
re-adjudicate it.

**Rule-7 scope sentence.** This campaign decided exactly the pre-registered
questions "C(U_mws59) <= 14" over the hash-pinned 389-element aux-1 universe at
T = 14, "C(V_mws59) <= 13" over the hash-pinned 366-element aux-1 universe at
T = 13, and "C(V_mws59) <= 14" over the complete construction-rule
two-auxiliary pair universe (79,728 pairs) at T = 14, for the `mws59` factor
blocks in their fixed orientation (sigma^0, all-monomial data triple) under the
three-stage linear SLP model with the transposition accounting; it did NOT
search the U-side two-auxiliary space at T = 15, the V-side three-auxiliary
space, any other orientation or data triple of `mws59`, any other
decomposition, non-ternary alphabets, or other models — so **whether the true
minimum is 58 or 59 remains open**, and no universal claim about `mws59` beyond
this orientation follows.

**Still open.** `mws59` fixed orientation 58-vs-59 (one auxiliary slot beyond
the registered recursion: U aux-2 at T = 15 or V aux-3 at T = 15, plus a
synthesis if admitted); `stapleton60` [58, 60]; `paper55` sigma^1/sigma^2
[55, 57]; `laderman23_source_lock`; and the frozen session-3
transposition-check artifact flagged above.

## Session 11 (2026-09-02, agent `Mm3NextGap`): `stapleton60`'s gap CLOSED at 60 — the published scheme is optimal for its own fixed orientation; the frozen LB 58 was not tight; verdict FROZEN-NEGATIVE

Frozen campaign:
`campaigns/20260902T012056Z_25e1634d_f146ed17f287/` (verdict `FROZEN-NEGATIVE`;
pre-registration committed BEFORE any feasibility compute at git `46c7618`,
copied byte-identically into the run dir with provenance binding).

**Headline.** `stapleton60` sigma^0 moves from [58, 60] to **LB = UB = 60**:
total = C(U) + C(V) + C(output) = 16 + 16 + 28, every term exact. Stapleton's
published 60-addition circuit is optimal for its own fixed orientation, and the
frozen certified lower bound 58 — still a valid bound — is NOT tight, because
both side floors' d+1 values (15, 15) are unattainable.

**The reduction.** C(output) = 28 exactly (frozen floor@13 on Wfac gives
output >= C(Wfac) + 14 >= 28; the printed output stage re-expanded in-run is
exactly 28 and reproduces all 9 C-entries over the 23 products), so the whole
gap is the two single-auxiliary questions C(U) <= 15? and C(V) <= 15? at
T = d+1 = 15. Closure argument (prereg section 3): a 15-gate circuit for 14 tau
classes has at most one non-tau gate value (two would leave 13 gates for 14
classes); duplicate-class gates are deletable at zero cost because operand
signs are free; a zero-auxiliary 15-gate circuit would strip to the refuted
14-gate floor schedule. Hence the pinned aux-1 universes are the whole space.

**Result — dual instrument, unanimous on all 1,252 instances.**
U side: 0 of 428 admitted at T = 15 => C(U) >= 16, and the printed 16-gate left
stage (re-expanded exactly, reproducing all 23 U rows) gives C(U) = 16.
V side: 0 of 426 admitted => C(V) = 16 likewise.
Wfac side: 8 of 398 admitted at T = 14; the first admission was replayed as an
explicit 14-gate circuit and verified by exact expansion => C(Wfac) = 14 exactly
(`witness_Wfac_14gates.json`). Universes hash-pinned before compute
(AU(U) 428 `9cb18bc7...`, AU(V) 426 `b0198869...`, AU(Wfac) 398 `29f85ab1...`)
and re-derived in-run. Cost: **110.7 CPU-seconds of a pre-registered 6.0 CPU-h
cap**, single thread, `nice -n 15`.

**Controls, both directions, all in-run, all first-pass.** ACCEPT: 729/729
Brent over Z in `fmpz` AND pure int; frozen landscape row reproduced
(14/14/13 | F,F,F | 1152/1344/460 | total_lb 58); printed-SLP recount by exact
expansion = 16/16/28 = 60 with every stage reproducing its rows; an
extension-positive control on the INDEPENDENT sun56-V target positive in three
instruments (DFS, kissat 4.0.4, CaDiCaL 3.0.1); transposition preconditions
audited (23/23 nonzero rows per side, rank(Wfac) = 9 exact over Q).
REJECT: one-gate-deleted and operand-perturbed plants both REJECTED by exact
expansion on the verified 16-gate carrier; a weaker-T negative (tau(Wfac) at
T = 12) certified fresh; d-count planter fired. **37 DRAT->LRAT certificates,
all accepted by BOTH pinned checkers**, replayed again post-run. Independent
post-run audit: **7/7 blocks pass**.

**No transposition-derived witness was used anywhere** — the session-10 D1
defect (the frozen session-3 `transpose_check.py` is self-circular) is honoured
by construction: the C(Wfac) = 14 upper bound is a synthesized explicit circuit
and the output-stage 28 is the printed stage re-expanded.

**Rule-7 scope sentence.** This campaign decided exactly "C(U) <= 15?" over the
hash-pinned 428-element aux-1 universe at T = 15, "C(V) <= 15?" over the
hash-pinned 426-element universe at T = 15, and "C(Wfac) <= 14?" over the
hash-pinned 398-element universe at T = 14, for the `stapleton60` factor blocks
in their fixed orientation (sigma^0, all-monomial data triple) under the
three-stage linear SLP model with transposition accounting; it did NOT search
other data triples or sigma classes of `stapleton60` (the frozen census gives
their lower bounds only), other decompositions, non-ternary alphabets,
GL(3,Q)/GL(3,Z) sandwiches, other counting conventions, or the open `mws59`
58-vs-59 residue, and it establishes no universal claim about rank-23 3x3
additive complexity. `paper55`'s 55 remains the record, untouched.

**Ladder status after sessions 9-11 (fixed orientations, all exact):**
`paper55` 55 = 55 · `sun56` 56 = 56 · `mws59` in {58, 59} (LB 58 certified) ·
`stapleton60` 60 = 60. Two published schemes (Sun's 56, Stapleton's 60) are now
certified optimal for their own orientations; `mws59`'s published 59 is within
one addition of its certified floor.

**Still open.** `mws59` 58-vs-59 (U aux-2 or V aux-3 at T = 15 plus synthesis);
`paper55` sigma^1/sigma^2 [55, 57]; `laderman23_source_lock`; owner
re-adjudication of the flagged session-3 `transpose_check.py` artifact.

## Session 12 (2026-09-02, agent `Mm3NextGap`): the `mws59` bit DECIDED — its fixed orientation is EXACTLY 58, with a verified 58-addition circuit; the published 59 is one above its own orientation's optimum; verdict FROZEN-CERTIFIED

Frozen campaign:
`campaigns/20260902T013730Z_49d2f937_a672bbec349b/` (verdict `FROZEN-CERTIFIED`;
pre-registration committed BEFORE any compute at git `f1668cd`, amendment 1 at
`b26a5dc` before the amended compute, prereg copied byte-identically into the
run dir with provenance binding).

**Headline.** total = C(U) + C(V) + C(output) = **14 + 15 + 29 = 58**, every term
exact, so `mws59`'s fixed orientation (sigma^0, all-monomial triple) has minimum
**exactly 58** (LB = UB = 58) — closing the residual bit left by session 10. The
published 59-addition scheme is **one addition above the optimum of its own
orientation**: its left stage spends 15 where 14 suffice. No published claim is
contradicted (59 additions do suffice), and **58 is not a record** — `paper55`'s
55 stands untouched.

**The decisive artifact.** `assembled_58_endtoend.json`: the assembled scheme
(session-10's frozen 14-gate left circuit + the paper's printed 15-addition
right stage + the paper's printed 29-addition output stage) was expanded
**symbolically over all 81 monomials A_i*B_j**, and every one of the 9 outputs
equals its matrix-product bilinear form `C[i][j] = sum_l A[3i+l]*B[3l+j]`
EXACTLY (9/9), with the total recounted from the explicit gate lists as 58.
This check is independent of the DFS/CNF searches, of the factor-block Brent
check, and of any transposition argument.

**How the two published stages entered as certified upper bounds.** The paper's
printed Table 2 (t0..t6, u0..u5, 23 products with inline signed operand
expressions, v0..v8, 9 output combinations) was transcribed verbatim from the
pinned layout artifact, recounted in the fixed model (**left 15 + right 15 +
output 29 = 59**, matching the paper and the frozen session-3 tally), and then
matched to the Table-3 factor blocks by an exact product correspondence:
coordinate maps (id, id, id), a bijection on the 23 products, per-product signs,
with **9 products carrying a negated value compensated in the output column**;
all 621 entry equalities hold exactly. Lower bounds come entirely from frozen
session-10 evidence: C(U) = 14 exact, C(V) >= 15 (aux-1 0/366 at T=13 AND
complete aux-2 0/79,728 at T=14, dual instrument), C(Wfac) = 15 exact,
output >= 29.

**Amendment 1 (disclosed).** The first Route-P attempt recounted correctly but
found NO correspondence, because the registered matcher forced eps*eps' = +1,
forbade a compensating output-column sign, and assumed a shared 3x3 flattening.
Amendment 1 widened the space to sigma_A, sigma_B, sigma_C in {id, T} plus
independent per-product signs, accepting only on exact equality of every entry,
and strengthened the wrong-correspondence plant to the widened space. The failed
narrow attempt is preserved verbatim and nothing from it is used. In the event
only the eps*eps' = -1 case was needed; the accepted coordinate maps are all
identity.

**Route S registered but never entered.** The fallback aux-3 existential at
T = 15 was registered with its closure argument and an enumeration-only
projection of **6.30e7 ordered chains at a measured 59.2-120 ms/instance =
173-2,100 CPU-hours**, i.e. two to three orders of magnitude beyond the 4.0
CPU-h cap; it was therefore registered as a prefix instrument with a
prefix-statement-only negative rule and a positive-only scouting tier. Route P
settled the bit, so no prefix was run, no sampling occurred, and no Route-S
claim exists.

**Controls and audit.** ACCEPT: 729/729 Brent in `fmpz` AND pure int; frozen
landscape row reproduced; BOTH frozen session-10 witnesses re-verified by exact
expansion; tri-instrument SAT control (DFS, kissat 4.0.4, CaDiCaL 3.0.1) on the
independent sun56-V target; transposition preconditions re-audited. REJECT:
one-gate-deleted plant, operand-perturbed plant, wrong-correspondence plant
(swapped labels in the widened space), frozen floor@12 re-certified INFEASIBLE
with a fresh dual-checker DRAT->LRAT certificate, and the d-count planter.
Independent post-run audit **6/6 blocks pass**. Cost: 2.2 CPU-seconds of a 4.0
CPU-h cap (the expensive work was session 10's 80,483-instance scans).

**Rule-7 scope sentence.** This campaign decided exactly the bit "is the mws59
fixed-orientation minimum 58 or 59?" for the sigma^0 all-monomial data triple
under the three-stage linear SLP model with transposition accounting, using
frozen session-10 lower bounds plus the verified printed Table-2 stages and an
end-to-end symbolic expansion of the assembled circuit; it did NOT search any
aux-3 space, other mws59 data triples or sigma classes, other decompositions,
non-ternary alphabets, GL(3,Q)/GL(3,Z) sandwiches, or other counting
conventions, and it makes no record claim.

**Exact fixed-orientation ladder after sessions 9-12** (published -> certified
exact minimum for that decomposition's own fixed orientation):
`paper55` 55 -> **55** (record, closed) · `sun56` 56 -> **56** (published
optimal) · `mws59` 59 -> **58** (published is +1) · `stapleton60` 60 -> **60**
(published optimal). Every rung of the published ladder now has an exact
per-orientation minimum.

**Still open.** `paper55` sigma^1/sigma^2 [55, 57]; `laderman23_source_lock`;
owner re-adjudication of the void session-3 `transpose_check.py` artifact; and
the global question of whether any rank-23 3x3 scheme beats 55 (untouched here).

## Session 13 (2026-09-02, agent `Mm3NextGap`): `paper55` sigma^1 and sigma^2 CLOSED at EXACTLY 55 — the record value holds at three orientations; verdict FROZEN-CERTIFIED

Frozen campaign:
`campaigns/20260902T015158Z_c0bea75d_761cb990fc40/` (verdict `FROZEN-CERTIFIED`;
pre-registration committed BEFORE any compute at git `bdf3a13`, copied
byte-identically into the run dir with provenance binding).

**Headline.** Both remaining `paper55` orientations collapse from [55, 57] to
**exactly 55**: sigma^1 = 14 + 14 + **27** = 55 and sigma^2 = 14 + 13 + **28**
= 55, each with an explicit circuit verified stage-by-stage by exact expansion
and END-TO-END by symbolic expansion over all 81 monomials (**9/9 outputs, zero
mismatches, both orientations**). Their "best known 57" rows are superseded by
2. This is the record VALUE at two further orientations of the same tensor —
**not a new record**: sigma^0's 55 remains the record and nothing sub-55 is
claimed.

**Why almost no search was needed.** Gate B's frozen results already make the
sigma^0 side costs exact (C(U) = 13, C(V) = 14, C(Wfac) = 14). The structural
step is the **proved free-relabeling symmetry**: a fixed permutation of the
input coordinates carries a circuit to a circuit of equal cost (rename free
input wires), so `C(T(X)) = C(X)` for the per-product transpose T — re-checked
numerically in-run (d and floor verdict of T(X) equal those of X for all three
blocks). With `sigma^1 = (V, T(W), T(U))` and `sigma^2 = (T(W), U, T(V))` from
the frozen `gatec_sweep.sigma_orbit`, every side cost is therefore exact and the
only open quantity per orientation was whether its output stage attains its
bound (27 and 28).

**What was constructed.** (i) A **14-gate Wfac circuit synthesized in-run** by
the aux-1 dual instrument over the hash-pinned universe |AU| = 391
(sha `1a3e66e9…`, two independent enumeration loops agreeing), admitted at
instance 307 and verified by exact expansion — so **no transposition-derived
Wfac witness is used**, honouring the session-10 D1 lesson. (ii) The printed
`LEFT_SLP` (13) and `RIGHT_SLP` (14) reused under free input renaming, each
re-verified. (iii) Each output stage built by **reverse-mode transposition** as
explicit adjoint accumulations with additions counted exactly, then checked
entry by entry (23x9 equalities) to reproduce its orientation's W-combination
map. (iv) End-to-end 81-monomial expansion per orientation.

**Controls and audit.** ACCEPT: 729/729 Brent for all three orientations in
`fmpz` AND pure int; printed-SLP recount 13/14/28 = 55 with all targets
reproduced; frozen landscape rows reproduced (d/state counts (12,33), (13,116),
(13,66) and their cyclic images); relabeling symmetry re-checked; transposition
preconditions audited per orientation; tri-instrument SAT control (DFS + kissat
4.0.4 + CaDiCaL 3.0.1). REJECT: gate-deleted and operand-perturbed plants both
rejected; the frozen Wfac floor@13 re-certified INFEASIBLE with a fresh
dual-checker DRAT->LRAT (rc 0/0); d-count planter fired. Independent post-run
audit **6/6**, including a from-scratch re-derivation of BOTH assemblies (55
each, 9/9 outputs). Cost ~81 CPU-seconds of a 2.0 CPU-h cap.

**Two instrument defects, both caught by pre-registered checks (disclosed).**
D1: the first assembly relabeled circuits by rewriting stored values, and the
exact-expansion verifier correctly REFUSED the stage; fixed by threading the
input permutation into the evaluator. D2: with a permuted input assignment the
transposed stage's 9 outputs are indexed by input WIRES and so emerged in
permuted coordinate order — the pre-registered output-map equality check caught
it while the counts were already correct (27/28, totals 55/55); fixed by
reporting the transposed outputs in coordinate order (a free relabeling). The
intermediate INCONCLUSIVE verdict is preserved verbatim; no claim rests on
either attempt.

**Rule-7 scope sentence.** This campaign decided exactly the minima of the
`paper55` sigma^1 and sigma^2 orientations under the three-stage linear SLP
model with transposition accounting, using the frozen gate-B exact side costs
transported by the proved relabeling symmetry plus in-run synthesis,
construction and verification of the missing stage circuits; it did NOT search
any multi-auxiliary space, other decompositions, other data triples or
sandwiches, non-ternary alphabets, GL(3,Q)/GL(3,Z), or other counting
conventions, and it makes no sub-55 claim.

**Exact fixed-orientation ladder after sessions 9-13** (published -> certified
exact minimum): `paper55` sigma^0 55 -> **55** (record) · `paper55` sigma^1
(57 best known) -> **55** · `paper55` sigma^2 (57) -> **55** · `sun56` 56 ->
**56** (published optimal) · `mws59` 59 -> **58** (published +1) ·
`stapleton60` 60 -> **60** (published optimal). Every named orientation with a
frozen LB/UB gap in this target is now closed.

**Still open.** `laderman23_source_lock`; owner re-adjudication of the void
session-3 `transpose_check.py` artifact; and the global question of any rank-23
3x3 scheme below 55, which no campaign here has touched.

## Session 14 (2026-09-02, agent `Mm3NextGap`): `laderman23_source_lock` Phase 1 CLOSED — Laderman's own 1976 paper pinned, transcribed verbatim, and machine-verified 729/729; verdict FROZEN-CERTIFIED

Frozen campaign:
`campaigns/20260902T022208Z_8c279ad1_d0276dfee63d/` (verdict
`FROZEN-CERTIFIED`; pre-registration committed BEFORE any target compute at git
`425c8a5`, amendment 1 at `cb319be`, both copied into the run dir with
provenance binding).

**The long-open `laderman23_source_lock` item is resolved.** Sessions 8-13 all
recorded it as blocked with "no verified factors on disk". Phase 1 obtained the
**primary** source — J. D. Laderman, *"A noncommutative algorithm for
multiplying 3x3 matrices using 23 multiplications"*, **Bull. Amer. Math. Soc.
82(1), January 1976, 126-128**, the open-access publisher PDF, 236,490 bytes,
**sha256 `a1deb200f2e270b13b870bbcf8ce8c669d7a09a4acdeddc012df5d6c90042a5b`** —
pinned it in the run dir, transcribed all 23 products and all 9 output equations
verbatim from the page, and verified the paper's OWN stated criterion (p.128,
"729 nonlinear algebraic equations involving 621 unknowns"):
**729/729 Brent identities hold exactly**, in pure Python `int` AND
`python-flint fmpz`, the two arithmetics agreeing on every single identity.
Triple sha256 `522ba07f4f1784ac…`, alphabet exactly {-1,0,1}, 23 products, no
zero rows; stored as `laderman23_factors.json` in the frozen run dir.

**Provenance labels, both carried and not interchangeable.** That these formulas
are Laderman's published scheme is a **CITED-DEPENDENCY** on the pinned
open-access PDF (authorship is cited, never machine-verified). That the
transcribed triple is a valid rank-23 decomposition of the 3x3 matmul tensor
over Z is **MACHINE-VERIFIED** in-run. That the transcription matches the
printed page is **MACHINE-VERIFIED** twice: the 7 single-product terms and all 9
output supports equal the printed equations, and every line equals the committed
amendment's recorded read term-for-term and sign-for-sign.

**Controls, both directions, all passed.** A1 729/729 dual arithmetic; A2 the
frozen `paper55` triple re-verified 729/729 through the SAME code path (a
known-true object must be accepted); A3 printed-page fingerprint; A4
alphabet/rank audit; A5 (added by amendment 1) zero drift against the recorded
primary read. REJECT: R1 one U sign flipped (1 identity broken), R2 one product
deleted (7 broken), R3 two V rows swapped (35 broken), R4 planted count
assertion fired. Independent post-run audit: a from-scratch Brent loop plus an
independent **symbolic** expansion showing each printed output equation expands
to exactly `sum_k a_ik b_kj`. Cost: under 1 CPU-second of a 0.5 CPU-h cap.

**Reproduction cross-check (non-gating, corroborating).** An independent pinned
reproduction (Courtois-Bard-Hulme, arXiv:1108.2830 §2.4, sha256 `d768596b…`)
transcribes to a different-but-equivalent presentation: it also satisfies
**729/729**, and its nine output supports are **identical 9/9** to the primary's
printed equations while per-product signs differ (negating a product's operand
together with its output coefficients is free). Support agreement across two
independently typeset sources corroborates that the locked object is Laderman's
scheme and not some other rank-23 decomposition.

**Two defects, both caught by pre-registered controls, no claim affected.**
D1: the first run aborted at A1 with 727/729 (both instruments agreeing on 2
failures) and produced NO verdict — the `pdftotext` layer truncates the final
term of `m_3` and `m_11` at the page edge and I typed `b_33` for both, where the
scan prints `b_33` for `m_3` and `b_32` for `m_11`. Amendment 1 (committed
BEFORE the amended run) demoted the text layer, made the 300/400 dpi
rasterizations authoritative, recorded the full line-by-line read, corrected
that one coefficient, and added gating control A5. **The prereg's
FROZEN-NEGATIVE branch was deliberately not invoked: my typing failed, the
source did not.** D2: A5's own parser initially tokenized a trailing prose
annotation as formula terms and aborted; parser fixed, recorded read unchanged.
The aborted verdict and pre-amendment runner are preserved verbatim.

**Rule-7 scope sentence.** Phase 1 produced exactly one artifact — a
hash-pinned, Brent-verified factor triple with its provenance record. It claims
**no** addition count, stage cost, floor, optimality, or comparison to any
published operation count, and it does not verify authorship. Phase 2 (the
ladder) requires its own fresh pre-registration.

## Session 15 (2026-09-02, agent `Mm3NextGap`): `laderman23` Phase 2 CLOSED — the exact fixed-orientation minimum of Laderman's 1976 scheme is EXACTLY 62, against a printed 98; verdict FROZEN-CERTIFIED

Frozen campaign:
`campaigns/20260902T023856Z_28ff8823_43330ec64df9/` (verdict
`FROZEN-CERTIFIED`; pre-registration committed BEFORE any feasibility or
synthesis compute at git `7e3321a`, amendment 1 at `5e40fc6`, both copied into
the run dir with provenance binding).

**Headline.** For the factor triple locked and Brent-verified by Phase 1,
at Laderman's own printed (sigma^0, all-monomial) orientation:

    C(U) = C(V) = C(WFac) = 16 EXACTLY,  output stage = 30 EXACTLY,
    total = 16 + 16 + 30 = 62 EXACTLY   (LB = UB = 62)

**Laderman's printed basic form recounts to exactly `28 + 28 + 42 = 98`
additions** (recounted twice independently), so the certified minimum sits
**36 additions below the printed count**. This is *not* a correction of
Laderman: his paper states outright that "the number of additions in this
algorithm could be greatly reduced, but it is being given in its more basic
form". The campaign **quantifies that advertised slack exactly** — the
reducible amount is exactly 36, and 62 is the floor.

**The decided question.** `total = C(U) + C(V) + C(WFac) + 14` with the
transposition gap (preconditions audited: 23/23 nonzero rows per side,
`rank(WFac) = 9`). In-run `d = 14` on all three sides with **all floors
impossible** (2401 states each) => each side `>= 15`, total `>= 59`. The open
question was whether 15 is attainable, i.e. the aux-1 existential at `T = 15`.
By the closure lemma (at `T = 15` exactly one gate is non-needed, so the
auxiliary is a signed pair-sum of inputs ∪ needed classes) the universe is
finite: **hash-pinned at `|AU| = 408` per side BEFORE any feasibility compute**
(`a1de3df4`, `d3f38942`, `13be8045`), re-derived in-run by two independent
enumeration loops.

**Decisive evidence.** The **complete census of 3 x 408 = 1,224
dual-instrument instances admitted ZERO** auxiliaries — memoized subset-DFS and
a fresh slot-availability CNF (kissat 4.0.4) agreeing on **every** instance,
index coverage verified complete 0..407 per side — so every side is `>= 16` and
the total `>= 62`. **All 1,224 infeasibility results carry DRAT->LRAT
certificates accepted by both pinned checkers (1,224/1,224 rc 0/0)**; proof
files are retained for 16 representatives with every certificate's hashes and
checker codes recorded per instance (retention policy disclosed in the report).
For the upper bound, randomized greedy CSE (seed `20260902`, 240 restarts/side)
proposed circuits of which **240/240 per side passed exact-expansion
verification, 0 discarded**, the best being **16 gates on every side**; the
output stage was then **constructed** by reverse-mode transposition of the
verified `WFac` circuit at exactly **30** additions and checked against the
`WFac` map in all **207** entries; and the assembled scheme was expanded over
all **81 monomials** — **9/9 outputs equal the 3x3 matrix product, zero
mismatches**.

**Controls.** ACCEPT: A1 frozen triple (sha `522ba07f`, Brent 729/729 int+fmpz);
**A2 a known-true circuit accepted at its exact addition count** (the frozen
`paper55` printed circuits verified at exactly 13 and 14 gates, recount
`13+14+28 = 55`); A3 frozen landscape row reproduced; A4 relabeling invariance
(8 permutations/side, seed 17); A5 preconditions; A6 universe pins via two
independent loops. REJECT: R1 gate-deleted plant, R2 operand-perturbed plant,
R3 floor `T = 14` re-certified UNSAT with a fresh dual-checked certificate, R4
planted `d`-count assert, R5 (amended) a known-SAT creatable chain feasible in
**both** instruments. Cross-solver: CaDiCaL 3.0.1 agrees with kissat.
Independent post-run audit **8/8**, including an independent re-transposition to
30 additions and an independent end-to-end expansion. Cost **~7.5 CPU-minutes**
of a 2.0 CPU-h cap.

**Defect disclosed.** D1: control R5 was mis-specified in my own
pre-registration — I called "classes + 5 slack slots" trivially schedulable, but
this encoding's operand pool is `inputs ∪ tau`, so slack slots add no power and
the instance is UNSAT exactly because the floor is. **The encoder was right and
my control was wrong**; the first run aborted at R5 with no verdict, and
amendment 1 (committed before any amended compute) replaced it with a
construction-known-SAT creatable chain now required to pass in *both*
instruments. Aborted verdict and pre-amendment runner preserved verbatim.

**Route not run, quantified.** The aux-2 census at `T = 16` (~90,000
pairs/side, measured ~2.0 CPU-h/side) was never needed: the ladder closed one
level earlier because the aux-1 census gave `LB = 16` and the verified 16-gate
circuits met it. Recorded as a quantified untaken route; nothing is claimed from
it.

**Rule-7 scope sentence.** This campaign decided exactly the minimum of
Laderman's own factor triple at its printed orientation under the frozen
three-stage model with transposition accounting, over a completely enumerated
hash-pinned aux-1 universe. It says nothing about other orientations of
Laderman's tensor, other decompositions, non-ternary alphabets,
`GL(3,Q)`/`GL(3,Z)`, other actions or counting conventions, or the global
sub-55 question. Laderman's 62 is **above** the record 55 and is not a
competitor for it.

**Exact fixed-orientation ladder after sessions 9-15** (published -> certified
exact minimum): `paper55` sigma^0 55 -> **55** (record) · sigma^1 (57) ->
**55** · sigma^2 (57) -> **55** · `sun56` 56 -> **56** · `mws59` 59 -> **58** ·
`stapleton60` 60 -> **60** · **`laderman23` 98 (explicitly unoptimized) ->
62**.

**Still open.** Owner re-adjudication of the void session-3
`transpose_check.py` artifact; and the global question of any rank-23 3x3
scheme below 55, which no campaign here has touched.

## Session 16 (2026-09-03/04, agent `Mm3Sub55`): `sub55-newspace` CLOSED — the full ternary orientation orbit of `laderman23` certifies >= 59 with ZERO <=54-certifiable points; the named landscape is now complete over all five public rank-23 factor triples; verdict FROZEN-NEGATIVE

Frozen campaign:
`campaigns/20260904T015105Z_49ed8737_988c6df478f1/` (verdict
`FROZEN-NEGATIVE`; prereg `sub55_newspace_pre_statement.md` committed at
git `e25d23e` BEFORE any compute, byte-identical copy + provenance in the
run dir).

**The decided space.** Sessions 5-8 swept the full ternary-unimodular
orbit $(P,Q,R)\in T^3$ ($|T|=6960$) x sigma-orbit for four factor triples;
`laderman23`'s verified factors did not exist until session 14 and session
15 decided only sigma^0. This campaign enumerated and decided the last
unswept named slice: $3\times6960^3 = 1{,}011{,}460{,}608{,}000$ raw
scheme-instances, through the proven 145-node factorization (universe
hash-pinned `bed89ca17f5bd868...` from two independent enumeration shapes)
to 34,944 admitted data triples (11,648 per sigma class, equal counts and
equal LB histograms exactly as sigma-conjugation predicts), each decided.

**Result.** Minimum certified total LB = **59** in every sigma class,
unique at the all-monomial data triple (125,125,125) — reproducing the
frozen session-15 sigma^0 anchor — with the next tier at 61 (12 triples
per class), max 86, and **zero rows <= 54**. The +14 transposition gap's
activity precondition held on all 34,944 triples (zero gap-drops).

**Instruments (dual, upgraded over sessions 7/8 which were DFS-only).**
Every positive pair-side instance (2,736 per sigma; 8,208 total) decided by
complete memoized subset-DFS AND kissat 4.0.4 on the slot-availability CNF
at $T=d$, with **100% instrument agreement**; all 7,029 floor-impossibilities
carry DRAT->LRAT proofs accepted by BOTH pinned checkers (drat-trim +
lrat-check, frozen 2e3b2dc snapshots); CaDiCaL 3.0.1 concordant on
representatives. All `.lrat` retained; full per-instance hash ledger.

**Controls, both directions.** ACCEPT: Brent 729/729 int+fmpz on the
hash-pinned triple (sha `522ba07f...`); 12 frozen anchor rows reproduced;
known-true circuits accepted at exact counts (paper55 13/14/28 = 55;
session-15 witnesses 16/16/16, recount 62); relabeling invariance 30/30;
145^2 full-row + 400-random direct-vs-table controls. REJECT: gate-deleted
plant (coverage), operand-perturbed plant (gate value), planted UNSAT floor
at T=14 with fresh dual-checked certificate + known-SAT chain, planted
monomial anchors, wrong-action plants. No transposition-derived witness
used. Independent post-audit: `postcompute_audit.py` — 138 checks,
0 failures (fresh survivor re-enumeration, histogram identity across
sigma, ledger coverage exact, retained-proof re-verification, 90 fresh
dual-DFS pair re-decisions).

**Defects (pre-verdict, preserved).** Two launch failures (a path bug; an
A3 target-set error — witnesses must be checked against laderman23's own
locked rows, up-to-sign per the free-sign-change semantics) and one
instrument repair (kissat's positional DRAT argument had been dropped;
aborted per prereg, fixed, pre-flighted). Attempt logs retained
(`runner_log.attempt*.txt`, `nohup_runner.attempt*.out`).
Post-freeze owner audit corrected one display-only arithmetic typo:
$3\times6960^3=1{,}011{,}460{,}608{,}000$.  The frozen `report.md` prints
$3{,}034{,}381{,}824{,}000$ after multiplying by three twice; the
pre-registration, `campaign_result.json`, reduction identity
$3\times145^3\times48^3$, and every enumerated count use the correct value.

**Rule-7 scope sentence.** This campaign swept exactly `laderman23` under
the full ternary-unimodular product and the 3-element sigma orbit under the
proven automorphism family; it does not sweep any decomposition outside the
five named public ones, non-ternary alphabets, GL(3,Q)/GL(3,Z) beyond T,
actions outside the proven family, anti-cyclic BA swaps, or upper-bound
synthesis; it establishes **no universal no-54 claim**. An LB is a lower
bound only: no circuit or optimality claim is made.

**Where the frontier now stands.** The named landscape — five public
rank-23 factor triples x full ternary orientation orbits x sigma — is
certified >= 55 end to end, with zero <=54-certifiable points. Remaining:
(i) `laderman23-orbit-min` — the orbit's exact minimum (preregistered
next; expected 62 exact via aux-1 pushes on the 36 LB-61 triples);
(ii) public decompositions beyond the five named (need primary
source-locks); (iii) the genuinely global question, for which candidate
(c)-style rank-profile bounds were CLOSED here as provably non-decisive
(universal per-side bound sharp at rank-9 => total >= ru+rv-4 <= 42 < 54).

## Sixth public source locked; census capped — CLOSED (Main, 2026-09-04)

[Run 20260904T035648Z_a4bf0c37_d790a569af43](campaigns/20260904T035648Z_a4bf0c37_d790a569af43/VERDICT.md)
is **FROZEN-INCONCLUSIVE** for its full census, with an independently verified
sixth-source lock: **S1_smirnov_repo_139**, 23 products, 139 nonzeros,
ternary, sigma-zero non-axis projective-class counts **(20,14,14)**.
The fixed-order selection no longer depends on ternarity; ternarity only
chooses its registered census path.

A driver-free source replay passed all **58 checks**, including all 729
standard tensor coefficients, the classical positive/deleted-product negative
controls, and exact invariant separation from all five frozen triples under
the registered equivalence action. This is not arbitrary linear-isotopy
classification or a new addition-count record.

The sampled 64 MiB growth cap stopped the census with **520 pair decisions**
and **zero complete safe orientation rows**. Peak sampled host CPU was
41.98%, group RSS 127.8 MiB; actual growth was 64.133 MiB. The prefix is
preserved but has no independent full-census audit or negative promotion.
All **1,872 frozen checksums** pass. Prior failed attempts and stale smoke
artifacts were preserved, not deleted.

A successor needs a preregistered resource-compatible census/evidence
strategy and independent replay before reusing the saved prefix. The former
five-source results remain scoped to those five; this sixth source's full
orientation landscape is still unresolved. No universal no-54 claim follows.


<!-- cs-autonomy:20260905T043238-e8dce0efdffe:1 -->
## Autonomous mm3 closeout — 2026-09-05 05:01 UTC

**FROZEN-CERTIFIED**. REVIEW of H-MM3-RESUME-BOUNDARY-1: VERIFIED within its registered scope (resumable ledger mechanism + retained 520-pair prefix integrity; no census progress, no addition-count or no-54 claim). Independent verifier review_verify.py (own code, primary driver read only as data; sealed sha256 cc1b781e…, one receipt rc 0, 2.0 s wall, 8 KiB growth) reproduced every registered finite claim: R1 prefix sha256 0adc781e…0e09, 520 unique keys = sigma|side|a|b, U/V/W 256/224/40, 8 SAT / 512 UNSAT, all UNSAT rcs 0 with lb=d+1/floor false, SAT cert null with lb=d/floor true, all 512+512 present cnf/lrat files hash to recorded values (8+8 missing = the SAT records), 0 anomalies. R2 own chain implementation of the documented H_0/H_i gives H_520 = 40887cc0…134e (= outA = outC) and H_200 = 8b875d99…8aba (= ckpt_200.json), ckpt last_key/tallies equal my recomputation over lines 1..200, outA.final == outC.final, six primary tampers rejected under the documented rules. R3 arithmetic reproduced exactly: cnf 44,758,113 B + lrat 17,139,872 B = 119,035 B/decision vs 724 B/decision ledger; 560 vs 92,646 decisions per 64 MiB. R4 triple ternary, 139 nonzeros (49,45,45), all 729 Brent identities hold in exactly one orientation (roles 0,1,2; my convention needs no transpose, the primary's convention reports transpose — a labeling difference, not a contradiction). Beyond the primary: R5 own pure-Python LRAT checker verified 5/5 deterministically sampled UNSAT pairs (0|U|20|1, 0|U|118|33, 0|V|40|38, 0|V|107|122, 0|W|77|115) independently of the recorded return codes; each proof is a single RUP step, i.e. these pair instances are unit-propagation refutable. Two findings: (a) the primary's open note is resolved — the locked triple hash cf9d004c…332c IS reproduced by the frozen lock's canonical form json.dumps({"U":U,"V":V,"W":W}, sort_keys=True) (independent_source_lock_audit.py lines 294–297); successors should pin this form. (b) Design gap, not a claim failure: the primary's checkpoint validator checks tallies only by n and sum consistency, so a compensated tamper (sat+1, unsat-1) passes validation (compensated_tamper_caught_by_primary_rules=false); the prereg wording 'refuse on any tally mismatch' is only met for the six tested tampers. A successor must recompute tallies from ledger lines 1..n (the chain already binds those lines, so this is cheap). Evidence: review/evidence/review_out.json, receipt compute-01.

Evidence: [`mm3/campaigns/20260905T045503Z_bbba5a52_6c94aeba1afd`](campaigns/20260905T045503Z_bbba5a52_6c94aeba1afd/autonomy/finalize.json). The label applies only to the registered claim; no broader frontier improvement is implied.

Next registered-work proposal: Harden the resume boundary before any census successor: validator must recompute side/sat/unsat tallies from ledger lines 1..n (chain-bound) instead of sum checks, and pin the canonical triple serialization json.dumps({'U','V','W'}, sort_keys=True) that reproduces cf9d004c...332c; certify with the compensated (sat+1, unsat-1) tamper plus the six existing tampers on the retained 520-line ledger.


<!-- cs-autonomy:20260905T091934-6caca9f7f32f:1 -->
## Autonomous mm3 closeout — 2026-09-05 10:15 UTC

**FROZEN-CERTIFIED**. REVIEW of H-MM3-ENCODER-FIDELITY-1: VERIFIED within its registered scope (encoder fidelity on every retained hash oracle + sigma^0/1/2 instance-universe registration; no solving, no census progress, no addition-count or no-54 claim). Own pure-Python verifier review_encoder_replay.py (sha256 4f95c925…49c5, git c6de5673) written from the frozen sources read as data (new_decomp_offdiag 143-360/757-820/905-924, gate_b_floor 38-74, orbit_min_run 244-331/477-480/2270-2276/2330-2363/2506-2523, GAP=14) with different design (flat 9-tuples, adjugate inverse, orbit closure then sort, early-exit tables, index-merged rep lists); primary gate script neither imported nor copied, its registration JSON read only as data. One compute receipt rc 0, sealed_intact, 36.2 s charged, peak RSS 30 MB, growth 250 KB. All 17 preregistered booleans true: R1 census 19683/11808/6960/48/4656/4608, 145 orbits x 48 partition 6960, unique factorisation, mono 125, pin sha256 prefix bed89ca17f5bd868 = log; R2 tables 256/224/128, 224/128/256, 128/256/224 with 1024 survivors each (dense count agrees, lexicographic, mono-node counts 735/259/29/1 in all three classes, direct-sandwich sample + negative grid consistent, sigma^0/1 = log); R3 attempt-3 pair file = 660 records (608 sigma^0 then 52 sigma^1), the 608 keys equal my ordered sigma^0 pair-key list in FULL ORDER (closes the primary's 520-prefix-only order check), the 52 equal my first 52 sigma^1 keys, 0 mismatches over cnf_sha256/cert.cnf_sha256/vars/clauses/d/lb/cert-presence/agree; R4 attempt-4 = first 520 sigma^0 keys in order, 0 mismatches, decisions equal attempt-3; R5 1024 row keys = my sigma^0 survivor keys in order, side_lbs = the three pair lbs, total_lb = sum+14 for all rows (min 63; 0 rows <= 54 in retained data only); R6 primary registration reproduced from my own lists: pair_keys_sha256 929d5d46…bd29, row_keys_sha256 6d4367a7…69b8, combined 96c3ef1a…2488, per-class instance tables equal row by row (3x608) with instances_sha256 1af9f689…ed9c / 894f316d…cca3 / 9cddb483…8330, identical d histograms; the primary's two false booleans are exactly the two that assumed 608 records, so its 'wrong count prior' explanation is confirmed and its PASS-on-every-oracle conclusion is correct. Informational (no claim): exact CNF-hash overlap between sigma classes is small (|s0∩s1|=28, |s0∩s2|=16, |s1∩s2|=8 of 608), so the candidate 'sigma permutes the instance universe' idea is not supported at identical-CNF level and needs a class-set canonical form to test.

Evidence: [`mm3/campaigns/20260905T100124Z_4e1069d8_da4f5d1086e0`](campaigns/20260905T100124Z_4e1069d8_da4f5d1086e0/autonomy/finalize.json). The label applies only to the registered claim; no broader frontier improvement is implied.

Next registered-work proposal: Driver-free sigma^1 pair census chunk 1 (H-MM3-SIGMA1-PAIRS-1): sealed pure-Python encoder (now doubly validated, 1180 oracle hashes) + in-sandbox pysat at T=d for the 608 sigma^1 instances; SAT models verified by clause evaluation, UNSAT accepted only after an own UP/RUP refutation replay; resumable ledger of (key,d,cnf_sha256,sat,lb) in <=300 s / <=8 MiB chunks, with the 52 retained sigma^1 decisions (lines 609-660 of attempt-3) as cross-instrument controls that must agree.
