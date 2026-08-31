# Pre-statement — campaign `2026-08-31T08-26-54Z_hgateMechanism`

This file is committed as the FIRST file of the campaign, before any compute.
Everything below is pre-registered: hypothesis classes, probes, anchors,
adjudication thresholds, and the rule-7 scope. Nothing here may be edited after
compute begins; changes go to the campaign report as recorded retractions or
corrections, never by silent rewrite of this file.

## 0. Question under test

H-GATE (owner hypothesis, frozen README 2026-08-30): **at every instance run on
this target, the weight-≤3 window is non-empty iff min_i s_i ≤ 3** (the window
over point supports S of weight |S| ≤ 3; V = the RS tensor code at
t = (1,1,1)). 20 consistent rows so far (13 gate-B frozen + 7 landing P2
non-empty rows), 19 of them retrospective; the single prospective read is PN10.
The 3 empty rows (PN4, PN5, PN10) and the 17 non-empty rows are all
Λ = Id, t = (1,1,1).

This campaign attacks the MECHANISM, priority order:

1. **Derive/prove the weight-≤3 window-emptiness mechanism** as a theorem or a
   precisely characterized obstruction (what the dual product-code structure
   cannot deliver, and why);
2. only if the mechanism resists: the weight-≥4 support-window census at a
   min(s) ≤ 3 instance (~11,700 s at N=140, owner estimate; we have not
   re-priced it);
3. test the s_i = 1 corner, untested in both directions.

Per the resource contract (one low-priority CPU-bound heavy process at a time),
item 2 is HELD and will only be started after KgBandClose's announced heavy
pass completes; preregistered order below is preserved while waiting.

## 1. Setup and conventions (stated alongside every map, per repo rule 7-d)

- Inst(q, s, t=(1,1,1), lam=None): N = s0 s1 s2 points of S_0 × S_1 × S_2 ⊂
  (F_q^×)^3; flat index `idx = (a0 s1 + a1) s2 + a2` (a_i indexes S_i in
  ASCENDING element order; `src/rs.py` `Inst.idx`).
- A direction-i LINE fixes the other two coordinates; at t_i = 1 it carries
  the 1-dim space spanned by its all-ones vector. Concretely V (the code) is
  the F_q-span, over all directions i and all lines l of direction i, of the
  indicator 1_line (lifted flat to F_q^N): V = span{ 1_{line_i(x)} :
  i ∈ {0,1,2}, x = 0..N/s_i − 1 }.
- **Derived shape of V at t=(1,1,1), Λ=Id** (derivation to be machine-checked,
  see A1): every v ∈ V has a UNIQUE representation
  v(a0,a1,a2) = x0(a1,a2) + x1(a0,a2) + x2(a0,a1), and conversely every such
  triple-sum is in V. Equivalently V = ker D0D1D2 where D_i is the one-step
  discrete difference along axis i acting on F_q^N (D_i kills exactly the
  summands that do not depend on coordinate i). Candidate identity:
  **dim V = N − (s0−1)(s1−1)(s2−1)**, to be machine-checked on every instance
  touched in this campaign (anchor A2). Immediate cross-check that PREDATES
  any compute of this campaign: the frozen unit (5,(2,2,2),(1,1,1)) has
  5^7 − 1 = 78,124 nonzero V vectors, i.e. dim V = 7 = 8 − (2−1)^3 — the
  identity holds there. At (31,(2,3,5)) it predicts 30 − (2−1)(3−1)(5−1) =
  30 − 8 = 22. [Rule-5 inline retraction, pre-compute: this document's first
  draft in-head-anchored "V-dim 82" at (31,(2,3,5)) — wrong; N = 30, not 90.
  Corrected here BEFORE any compute or freeze; engine V_dim at A1 remains the
  machine verdict.]
- General Λ_i at t_i = 1: L_i(C_i) = { λ_i(a_i) 1(x0,x1) } pointwise-scaled.
  Divide every codeword pointwise by Λ0(a0) Λ1(a1) Λ2(a2) ≠ 0 (all λ_i entries
  are nonzero by Inst's assert, all S_i ⊂ F_q^×, so the product is nonzero
  and none of the summand denominators vanish): this maps V_Λ bijectively to
  V_Id, with x_i ← x_i ∏_{j≠i} Λ_j(a_j). **Supports are preserved by the
  pointwise division (nonzero scaling of each coordinate's value).**
  ⟹ Every emptiness/non-emptiness statement proven for Λ=Id holds verbatim
  for all Λ. THEOREM-level claim (F_q-linear bijection, pointwise nonzero
  rescaling); supports and weights are per-coordinate-value-invariant.

## 2. The ⇐ direction is a theorem; exact hypothesis class

**H-GATE-⇐ (theorem, one line):** if some s_i ≤ 3 (i.e. min(s) ≤ 3), then the
weight-≤3 window is non-empty. Proof: the indicator of a single direction-i
line is in V (it IS the line basis element, t_i=1), has weight s_i ≤ 3, and so
for S = that line's s_i points, V ∩ F^S ≠ {0}. Every such support yields ≥1
window witness. [delta = s_i for the single-line witness is NOT claimed
here — that would need an exhaustive cheaper-support sweep per witness; only
non-emptiness is asserted, matching H-GATE's semantics.]
Covers s_i = 1 too: a 1-point "line" is a point indicator, weight 1, trivially
in V (the line through a point has length 1 when s_i = 1).

Note this theorem's proof mechanism (single line indicator) is the
min(s)-controlled object; the 3/5 witnesses at (2,3,5) are NOT single lines
(single 3-line gives wt 3 = δ 3, ratio 1; the 3/5 rows have wt 3, δ 5,
a fundamentally different support). Do not conflate. The theorem only asserts
non-emptiness, matching H-GATE's semantics.

## 3. The ⇒ direction — rigidity, precisely what is open

**Claim under proof attempt (min-form):** if s_i ≥ 4 for all i, then
V ∩ F^S = {0} for every point support S with |S| ≤ 3.

An F_q-linear rigidity statement. Proof strategies pre-registered here:

- (a) Pigeonhole on a fiber: pick a coordinate (say a0-axis) fiber structure
  induced by S; if any fiber contains ≥2 points of S, restrict to the
  a0-lines through them. There are at most |S| ≤ 3 relevant fibers per axis;
  the x1, x2 summand is constant in a0 on all but ≤ 3 fibers. Second-difference
  pigeonhole fails at s_i = 4,5 (excluded sets grow past available points);
  needs a global argument or a genuinely false claim near singular shapes.
- (b) ker D0D1D2 small-support form: v ∈ ker D0D1D2 with wt(v) ≤ 3 ⟹ v = 0
  when all s_i ≥ 4. Equivalent formulation (affine-line moment vanishing):
  for every axis-0 line not fully supported, its line sum vanishes; and
  along-axis D_i v vanishes everywhere. v supported in ≤3 points: each of the
  ≤3 points' axis-mates must have certain moment cancellations; at s_i ≥ 4 the
  "over-determination" count exceeds the free coefficients available, which is
  where the theorem should live. [DERIVED plan, not a proof.]
- (c) Dual identity (dimension count): dim(V ∩ F^S) = dim V − |S^c| +
  dim(V^⊥ ∩ F^{S^c}) where V^⊥ = {h : every axis-parallel line-sum of h is 0}
  (orthogonal to V's generator set = the line indicators, since
  ⟨h, 1_line⟩ = line sum of h). Rigidity = dim(V ∩ F^S) = 0 for |S| ≤ 3, i.e.
  dim V − (N − 3) + dim(V^⊥ ∩ F^{S^c}) = 0, i.e.
  dim(V^⊥ ∩ F^{S^c}) = (s0−1)(s1−1)(s2−1) − 3, using dim V^⊥ =
  (s0−1)(s1−1)(s2−1) = N − dim V (the difference-operator row count).
  This is written here as a CHECKABLE identity to be anchored numerically
  BEFORE being used as a proof device; if the identity or its numeric anchor
  disagrees, the section is a retraction, recorded in the report.

## 4. Pre-registered probes (order fixed; executed in this order)

The ⇒ direction's hypothesis class has never been probed on REPEATED orders.
All frozen evidence has min(s) ∈ {2,3} (non-empty) or three DISTINCT s_i ≥ 4
(PN4: (4,5,7), PN5: (4,5,6), PN10: (5,7,8) — all pairwise distinct). If the
crude min-form were true, the following must ALL come out with EMPTY weight-≤3
windows. If ANY lands non-empty, the min-form of H-GATE is DEAD as stated.

| # | q  | s          | N   | # supports (w≤3) = C(N,1)+C(N,2)+C(N,3) | rationale |
|---|----|------------|-----|------------------------------------------|-----------|
| P1 | 17 | (4,4,4)   | 64  | 64 + 2,016 + 41,664 = **43,744**         | all three orders equal, 4\:16 ✓ |
| P2 | 41 | (4,4,5)   | 80  | 80 + 3,160 + 82,160 = **85,400**         | s0 = s1 repeated, both 4 |
| P3 | 61 | (4,5,5)   | 100 | 100 + 4,950 + 161,700 = **166,750**      | s1 = s2 repeated, both 5 |

s_i | q−1 feasibility: 4\:16 ✓; (4,4)|(40−?) — 4,4,5 pairwise-not-required
(subgroup orders need only divide q−1 individually: 4\:40 ✓, 5\:40 ✓);
4,5\:60 ✓. All three are pairwise-coprime-violating shapes IN PART (P2, P3) or
entirely (P1); but coprimality is known not to drive emptiness (PN5 non-coprime
is empty; (3,5,16) coprime is non-empty) — these probe REPEATED ORDER, not
coprimality.

## 5. The s_i = 1 corner (both directions untested by machine)

Instances pre-registered:

| # | q | s          | N | prediction from H-GATE-⇐ |
|---|---|-----------|---|---------------------------|
| C1 | 5 | (1,2,2)  | 4 | every support non-empty; min-ratio 1/2 (wt 2 line) or pt-indicator |
| C2 | 7 | (1,2,3)  | 6 | same |

These are shape probes outside the frozen-instance design conventions
(PN8/PN9-style repeats precedent); flagged as such. `subgroup(q,1,g)` returns
[1] (asserted s ≥ 1, loop never enters); line machinery handles it (line of
length 1; line_basis_flat is one vector). Machine-verify the ⇐ assertion
(new claim: s_i=1 is INSIDE the ⇐ arm, currently INFERENCE) and check whether
the ⇒ arm... vacuous here since min(s)=1 ≤ 3. The corner's real content: does
anything BREAK in the code at s_i = 1 (a robustness hazard for the frozen
instrument if ever reused downstream), and is the min-form's boundary exactly
"min ≤ 3" including min = 1.

## 6. Startup anchors (must pass before any probe may be reported)

A0. **Kernel re-validation (assignment mandate):** run the generalized census
kernel on frozen non-empty instances and reproduce
PN6 (43,(2,6,7)) → 3486 nonzero supports and PP2 (43,(2,3,4)) → 917 nonzero
supports, byte-equal to the frozen README numbers (384-386). Mismatch ⟹ HALT,
hub Main, report the discrepancy, do NOT proceed to probes.
A1. **V-structure check:** on (31,(2,3,5),(1,1,1)), verify the derived V
characterization: dim V == N − (s0−1)(s1−1)(s2−1) == 30 − 8 = 22, and on a
second instance, (13,(2,2,4),(1,1,1)): dim V == 16 − (2−1)(2−1)(4−1) =
16 − 3 = 13.
Also verify V ⊇ summands (dim of the naive 3-summand lift space) and
V ⊆ ker D0D1D2 by checking dim ker(D0D1D2) == dim V via an exact F_q rank of
the stacked difference operator. [Sub-second; part of the cheap controls.]
A2. **Support-count formula:** on every instance censused (probes and
anchors), assert at STARTUP, |all_S| == C(N,1) + C(N,2) + C(N,3) with the
exact integers from the table in §4 / frozen rows. Asserted before the loop
runs (repo rule 17b).
A3. **Instrument soundness counterfactual (controls give numbers meaning):**
before any emptiness verdict, plant KNOWN-NON-EMPTY cases and confirm the
instrument catches them: on a min(s) ≤ 3 instance (the PP2 instance suffices,
3486/917 nonzero supports expected... use PN6/PP2 numbers as the control),
the census must return nonzero counts (and by A0 already does). On the s_i=1
corners (C1, C2), the prediction is every-S non-empty; if C1/C2 return a
single zero-support, the instrument contradicts its own ⇐-theorem input —
halt, the instrument is broken.
A4. **Per-census anchors (rule 17b):** each probe run re-derives |all_S|
inside the script and asserts against the pre-registered count in §4 (43,744
/ 85,400 / 166,750) before enumerating.

## 7. Adjudication (pre-committed, no post-hoc)

- Probe window EMPTY ⟹ min-form's ⇒ arm under the enlarged hypothesis class
  (all s_i ≥ 4, including repeated and equal orders) survives these probes;
  evidence label: the census is exact F_q (rref mod q, python-flint / numpy
  integers only), so emptiness would be MACHINE-VERIFIED **within the
  candidate class and the window** (weight ≤3 supports; normalized reduced
  basis; zero nonzero V∩F^S found). H-GATE acquires 3 new consistent rows
  (now 23), all RETROSPECTIVE relative to the mechanism question — the
  mechanism itself remains OPEN unless section §3's proof closes (see §8).
- [RULE-5 INLINE RETRACTION, PRE-COMPUTE] The line this replaces originally
  read: "Probe window NON-EMPTY ⟹ H-GATE min-form REFUTED (a counterexample
  by complete enumeration, MACHINE-VERIFIED as a witness, exact F_q). ACTION:
  hub Main IMMEDIATELY…". It is FALSIFIED BEFORE ANY COMPUTE by the
  independently audited owner quotient proof (this audit and both of its
  pushbacks are recorded in the campaign report): at t=(1,1,1),
  V = ker(π0⊗π1⊗π2) with d(C_i) = s_i−t_i+1 = s_i, so for w < min(s) the
  quotient images of any support are independent, forcing V ∩ F^S = {0} —
  the ⇐/⇒ biconditional is now a THEOREM at t=(1,1,1), covering s_i=1,
  equal orders and all Λ (quotient argument is Λ-blind since d is
  Λ-preserving; the pointwise-division bijection gives the explicit Λ=Id
  route). A non-empty P1/P2/P3 window can now only mean an INSTRUMENT BUG
  (halting condition), not a refutation. P1–P3 remain run as
  pre-registered CONFIRMATORY censuses.
  Scope retraction, same source: the general-t predicate is
  min_i(s_i − t_i + 1) = min_i d(C_i), NOT min(s); the frozen q=13
  t-profiles agree with min(s) only because another axis still carries the
  small distance. H-GATE-as-stated (min(s) form) is a t=(1,1,1)-slice
  THEOREM, valid exactly when each C_i is the stated (Λ-diagonal) GRS code;
  pairwise coprimality is irrelevant to this distance theorem.
- C1/C2 any-zero ⟹ contradicts the ⇐ theorem that must hold; the instrument
  is BROKEN — halt, debug, re-run A0-A4 before any further verdict.

## 8. What this campaign does NOT claim (rule 7 scope, stated now)

Window means weight-≤3 point supports S, candidates = normalized reduced
basis of V ∩ F^S (nonzero). Swept domain: exactly the 3 probe instances P1-P3
and 2 corner instances C1-C2, plus the A0 anchor instances PN6, PP2,
(31,(2,3,5)), (13,(2,2,4)), all at t=(1,1,1), Λ=Id, complete weight-≤3
support enumeration per instance. NOT searched: any other (q,s) triple;
supports of weight ≥4; non-basis field values inside dim>1 intersections
(candidates are restricted to the basis class); other t; other η; non-identity
Λ (but see §1's Λ-uniformity: at t=(1,1,1) the Λ=Id emptiness on these probes
extends to all Λ by derivation — CHECKED pointwise-division bijection, theorem
strength, not swept by machine); the global ρ_inst everywhere (nothing here
reports ρ_inst; every landing "ratio" is a rho^window bound statement, and
empty windows only certify emptiness, no ratio at all). "No improving direction
found among those tried" is not "none exists" — the scope sentence in the
report will name exactly this.

## 9. Prospective/retrospective labeling policy

Every structural claim from this campaign is labeled PROSPECTIVE if its
prediction was written into THIS file before the compute that tested it (the
P1-P3, C1-C2 outcomes, A1's dim-V identity — all pre-registered here), and
RETROSPECTIVE otherwise (any pattern noticed after results landed). The two
"observes" in §1-§2 (⇐ theorem, Λ-uniformity) are derivations made before any
probe compute; they are labeled DERIVED, machine-anchor pending.

## 10. Compute budget and heavy-slot contract

The three probes are priced ~11,700 s total at N=140 (owner estimate from the
frozen P2 pass; we have NOT re-priced — the actual 2026 price comes from
P1-P3's own measured rates). Per the repo resource contract (one heavy
CPU-bound process at a time, Omega owns the slot until KgBandClose takes it),
probes are HELD until that handoff; this pre-statement, the anchors A0-A4,
and the s_i=1 corners (tiny, sub-second each) are cheap exact controls and
run under `nice -n 10` with thread pinning. Order preserved: A0-A4 first,
then C1/C2, then TP1-3. Instruction is explicit: no reorder of the probe
priority is possible anyway (P1 < P2 < P3 by N).
