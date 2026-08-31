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
