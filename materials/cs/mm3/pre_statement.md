# Pre-statement — mm3 campaign 1 (gates A, B, C)

**Committed 2026-08-29, BEFORE the first computational run of gates A/B/C.**
Agent: mm3 (Mm3). Folder: cs/mm3/.

## 0. Primary reads performed before this statement

Fetched first-hand (page + full HTML where available), 2026-08-29:

| ID | Exact title (quoted from fetched page) | Role |
|---|---|---|
| arXiv:2607.28676 | "55 Additions Suffice for 3×3 Matrix Multiplication at Rank 23" (Karunaratne, Idamekorala) | main paper: SLP, tensor factors, optimality claim |
| arXiv:2508.03857 | "A 60-Addition, Rank-23 Scheme for Exact 3x3 Matrix Multiplication" (Stapleton) | ladder step 60 |
| arXiv:2601.05272 | "A Rank 23 Algorithm for Multiplying 3 x 3 Matrices with an Arithmetic Complexity of 59" (Mårtensson, Stankovski Wagner, Stapleton) | ladder step 59 |
| arXiv:2512.21980 | "A 58-Addition, Rank-23 Scheme for General 3×3 Matrix Multiplication" (Perminov) | ladder step 58; source of cr58_cn122 |
| arXiv:2604.27645 | "An Exact 56-Addition, Rank-23 Scheme for General 3×3 Matrix Multiplication" (Sun) | ladder step 56 |
| Perminov repo, commit 98ba522db92b74f1f8c561a78038ff3091356d73 | file `schemes/results/addition_reduced_ZT/3x3x3_m23_cr58_cn122_ZT_reduced.json` | public cr58_cn122 realization (fetched; stored in campaigns archive) |

All five arXiv IDs resolve; titles quoted verbatim from the fetched pages.
The README's "four predecessor record papers" all exist as listed.

## 1. What each record optimized, and under which convention

All five papers state their count in the **standard bilinear straight-line
model**: inputs free; a gate computing x+y or x−y costs 1; copies and sign
changes cost 0; the 23 products counted separately; total = left + right +
output additions.

- 60 (Stapleton): "reduce the additive cost ... from the previous records of
  61 (Schwartz-Vaknin, 2023) and 62 (Martensson-Wagner, 2025) to 60 **without a
  change of basis**". Ternary factors [Sun's Table 1]. No optimality claim.
- 59 (M-W-S): combines Mårtensson-Wagner optimization with Stapleton's, "still
  without a basis change". Ternary factors. No optimality claim.
- 58 (Perminov): automated ternary-restricted flip-graph search + greedy
  intersection reduction (CSE); total scalar ops 83→81. Ternary. No
  optimality claim. Public artifact = the JSON above, claimed additive cost
  58 in its normalized structure (right side shares 14 additions with left).
- 56 (Sun): cyclic-automorphism reorientation (U,V,W)→(V,W^T,U^T) of the 58
  scheme, then re-synthesis: 13+13+30 = 56. Explicitly disclaims optimality
  ("we do not claim that 56 additions is globally optimal"). Ternary, 729
  Brent identities over Z. Its output count uses "t atoms → t−1" description,
  equivalent to the SLP-gate model.
- 55 (Karunaratne–Idamekorala): factor-level exact synthesis of the SAME
  cr58_cn122 orientation: U:9→23 cost 13, V:9→23 cost 14, W transposed →
  output 28; total 55. **Claim is explicitly narrow**: "For the fixed oriented
  tensor displayed in this paper, the minimum sum of the three linear-stage
  costs ... is 55", via d(F) direction-counting + exhaustive search at the
  d(F) floor: C(U)=13 (floor d=12 impossible), C(V)=14 (floor 13 impossible),
  and W-side: any W^T circuit of L gates transposes to a W circuit of L−14
  gates, so L ≥ C(W)+14 = 28 since C(W)=14 (floor 13 claimed impossible).

Convention caveat recorded: compare counts only within this model; the 58
artifact JSON's cost decomposition (their "reduced": 58) is a normalized
convention, not the paper's left/right/output split.

## 2. Gate A — confirmation harness (729 Brent identities, exact Z)

- **Quantity**: for all i,i',j,j',k,k' ∈ {0,1,2}:
  Σ_{r=0..22} U[r,(i,k)] V[r,(k',j)] W[r,(i',j')] = 1 iff i=i', j=j', k=k'
  else 0, evaluated in EXACT integer arithmetic (`python-flint fmpz`).
- **Inputs**: the tensor factors (U,V,W) derived two ways and required to
  agree: (a) from the paper's straight-line program (Sections 4.1–4.3, 13 U
  gates, 14 V gates, 28 W gates + free aliases); (b) from the paper's printed
  factor arrays (Section 5.2, three 9-row blocks of 23 columns). The SLP
  counts are also recomputed: left additions = 13, right = 14, output = 28.
- **Pass**: all 729 identities exact, counts 13/14/28/55, SLP-derived factors
  ≡ printed factors. Label [REPRODUCED] (same published object, independent
  transcription [DERIVED from paper text]).
- **Fail/refutation**: any identity fails, or counts mismatch ⇒ STOP, message
  Main via hub with the failing instance BEFORE writing anything.

## 3. Gate B — decision encoding for the per-orientation optimality claim

Fixed input: the paper's tensor factors (frozen in section 5.2).

**Question B(U)**: min # addition gates in a linear circuit over {-1,0,1}-less
gates (gates x+y, x−y of previously available quantities; inputs =
A_0..A_8 free; copies/sign-flips free) computing the 23 left factors U_r.
Is it ≤ 12? Paper claims 13 (minimum). Analogous B(V) with ≤ 13 (paper: 14).
Analogous B(Wᵀ) with ≤ 27 for the OUTPUT map (paper: 28), i.e. the transpose
9→23 factor circuit of the 23→9 output map; the W claim is C(W)=14 +
transposition bound, so a 27 gate W^T circuit would also break the W claim.

**Model (transcribed exactly from the paper's model, Section 2)**: linear SLP
with 9 free input wires; each addition gate g = x ± y where x,y are previous
wires or their negations (free sign absorption); sign changes and copies are
free (so a wire may be used negated at no cost — equivalently gates work on
signed quantities). The circuit computes the 23 target vectors exactly.

**Encoding 1 — ILP (HiGHS)**: variables x_{d,g} ∈ {0,1} = "gate g outputs
value-vector d ∈ Z^9" over the finite set D of attainable integer vectors
with ||d||_1 ≤ K (K = 12 for U, 13 for V/Wᵀ; gates with ||d||_1 > K are
pruned since every gate is a signed sum of two earlier vectors and target
reachability is checked exactly), plus solver-attested upper bounds via
iterative refinement if ever needed. Constraint: every target vector t (23
of them) has a computed wire t = x_{t,g} (targets may also be inputs);
every used gate wire has its inputs strictly earlier; "used" gate ⇒ its
value is a signed sum of two previous used values ⇒ for each candidate
value d and gate g: x_{d,g} ≤ Σ_{d1,d2 signed sum to d, d1<d positions}
... (standard adder-chain constraints, all-integer coefficients). Linear
relaxation dual extracted for the LP lower-bound certificate.

**Encoding 2 — SAT/PB (python-sat)**: cardinality/PB constraints with an
explicit Tseitin-free registration path per value; totalizer or sequential
counter for "# gates used ≤ T" with T = claimed−1 (12 / 13 / 27). Any SAT
answer (a feasible circuit) is re-verified OUTSIDE the solver by direct
SLP evaluation (the gate list is the certificate and gate A's harness
checks it); an UNSAT answer is only a verdict when accompanied by a
certificate.

**Certificate policy (fixed in advance)**: UNSAT ⇒ must be certified.
Attempted chain, in order of preference: (i) VeriPB proof from python-sat
PB encoding (tool `pysat.pb`–exportable .opb → whether a proof is produced
by our solver stack; check availability of VeriPB verifier, pin version);
(ii) if no proof-capable PB solver is available, the ILP route: exact-rational
LP relaxation dual (HiGHS dual certificate, rationalized + independently
re-multiplied against the constraint matrix), which for an infeasible
integer program still bounds the integer optimum: if LP optimum (#gates) >
T with a valid dual, then no integer circuit ≤ T exists **for the encoded
model**. Model-soundness caveat is recorded verbatim in the campaign README:
both encodings bound the SLP model; the certificate proves the MODELED
question. If encodings disagree with the paper's claim, MAIN MUST BE
NOTIFIED before writing a refutation anywhere.

**Verdicts**: feasible ≤ T ⇒ explicit gate list, re-checked through gate A
harness (as a factor circuit plugged into the 729 test), escalate to Main.
UNSAT at T with valid dual/proof ⇒ paper's per-orientation claim becomes
MACHINE-VERIFIED for that factor map at model level, label
MACHINE-VERIFIED, with the caveat "for the exact finite-gate-decision
encoding, not the free-algebra paper model; the paper's exhaustive-search
claim additionally covers unscheduled circuits, ours covers scheduled ones"
— and also escalate since a confirmation of optimality is worth recording.

## 4. Gate C — record attack / orientations

- **Fixed before launch**: sweep set = all orientations of cr58_cn122 that
  preserve rank-23 ternary structure, enumerated EXACTLY as:
  { row+column permutations of U and V (per row-major 3×3 structure),
    column permutations of the 23 product slots,
    transpositions (V,W^T,U^T) cyclic orbit } generalized as follows:
  we sweep the 6 orderings of "which factor block sits as left / right /
  output-transposed" × { row/column-permutation sweep NOT EXHAUSTIVE over
  the full S_9×S_9×S_23 group — budgeted }.
- **Budgeted sweep**: 6 factor-block rotations × per-rotation greedy CSE
  synthesis. NOT claimed exhaustive; report set exactly.
- **Per-instance solver timeout**: 900 s CPU per ILP instance, 600 s per SAT
  instance (nice -n 10, 1 thread). Total campaign wall-clock for gate C ≤ 8 h.
- **Decision**: total additions over all three stages treated as pure integer
  program; objective = sum of gate-usage indicator variables. Any total ≤ 54
  ⇒ new record, ship SLP + gate A re-run + own-orientation ILP, escalate to
  Main before claiming.
- **Solver versions pinned and recorded in the campaign manifest.**

## 5. Environment (fixed now, recorded in campaigns)

- Python: /Users/jinleic/jinleic-workspace/cs/.venv/bin/python 3.14.3
- python-flint 0.9.0 | numpy 2.5.2 | python-sat 1.9.dev15 | highspy (>=1.11)
 | z3 5.1.0 | sympy 1.14.0 | mpmath 1.3.0
- Resource policy: nice -n 10; single-thread pinning for BLAS/OMP via env
  vars written in run scripts.
- No additional packages are planned; any addition is recorded here + in the
  campaign manifest.

## 6. What can falsify this statement

- A gate A failure ⇒ published SLP/factors wrong ⇒ STOP + Main.
- A gate B feasible instance below the claimed floor ⇒ per-orientation
  optimality claim of the paper is FALSE for that factor orientation;
  certificate = explicit circuit list verified by gate A harness.
- A gate C total ≤ 54 ⇒ current published record is beaten; ship + escalate.
- Any other outcome (all floors confirmed) upgrades the claim to
  machine-checked optimality, one factor map at a time; report as such.

---

## Addendum (2026-08-30, agent `ProofLog`) — certificate route executed: DRAT/LRAT

Certificate policy section 3 anticipated "DRAT/VeriPB": the executed route is the first
arm of alternative (i)'s spirit — **DRAT proof logs from kissat 4.0.4, verified by
drat-trim (commit 2e3b2dc) and lrat-check** — all three floor questions certified.
Campaign: `campaigns/2026-08-30T031544Z_e0f3f117_c9df97a8bf3a/`. The SAT/PB and VeriPB
routes were NOT needed; RoundingSat/VeriPB remain untouched. Encoding-adequacy control
executed and passed as recorded in the campaign manifest (aux-1 control at d+1 = SAT,
paper witnesses admitted; raw floor control at d+1 = UNSAT because witnesses need one
aux intermediate — semantics of the floor question, not a bug).
