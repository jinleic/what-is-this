# Campaign report — gate `stapleton60-gap-total60` (run 20260902T012056Z_25e1634d_f146ed17f287)

## Question (pre-registered: pre_statement.md, sha c0203bb0…, commit 46c7618)

Decide the fixed-orientation minimum of `stapleton60` (Stapleton,
arXiv:2508.03857) inside the frozen gap **[58, 60]** — certified LB 58 from
the frozen landscape, published circuit 60.

## Answer — verdict FROZEN-NEGATIVE: the gap is CLOSED at 60

**LB = UB = 60 exactly. Stapleton's published 60-addition scheme is optimal
for its own fixed orientation, and the frozen certified LB 58 was valid but
NOT tight.** No 58- or 59-addition circuit exists for this orientation in the
three-stage linear SLP model.

    total = C(U) + C(V) + C(output) = 16 + 16 + 28 = 60

| stage | frozen floor | decision this run | exact value |
|---|---|---|---|
| C(U) | d = 14, floor@14 impossible (states 1152) | **0 of 428** single-aux extensions admitted at T = 15 (dual instrument, complete) ⇒ C(U) >= 16; printed 16-gate left stage re-expanded exactly | **16** |
| C(V) | d = 14, floor@14 impossible (states 1344) | **0 of 426** admitted at T = 15 ⇒ C(V) >= 16; printed 16-gate right stage re-expanded exactly | **16** |
| C(Wfac) | d = 13, floor@13 impossible (states 460) | **8 of 398** admitted at T = 14; first admission replayed as an explicit 14-gate circuit and verified by exact expansion | **14** |
| C(output) | >= C(Wfac) + 14 = 28 (transposition, preconditions audited) | printed output stage re-expanded: 9 v-gates + 19 combination additions = 28, reproducing all 9 C-entries over the 23 products exactly | **28** |

The closure argument (pre_statement §3) makes the two single-auxiliary scans
the WHOLE remaining question: a 15-gate circuit for a 14-class target set has
at most one non-tau gate value (two would leave 13 gates for 14 classes), a
duplicate-class gate is deletable at zero cost (operand signs are free), and a
zero-auxiliary 15-gate circuit would strip to a 14-gate floor schedule, which
the frozen floor decision refutes. So C(U) <= 15 ⟺ some a ∈ AU(U) makes
tau(U) ∪ {a} schedulable in 15 slots — decided negative over the complete
hash-pinned universe, and likewise for V.

## Search space, budget, measured cost

| side | d | \|AU\| (pinned) | T | instances decided | admissions |
|---|---|---|---|---|---|
| U | 14 | 428 (`9cb18bc7…`) | 15 | 428 | 0 |
| V | 14 | 426 (`b0198869…`) | 15 | 426 | 0 |
| Wfac | 13 | 398 (`29f85ab1…`) | 14 | 398 | 8 |

- Universes hash-pinned in the pre-statement (enumeration-only, before the
  prereg commit; **zero feasibility instances were decided pre-commit**) and
  re-derived + asserted in-run before each scan.
- **Dual instrument on every one of the 1,252 instances** (frozen complete
  memoized subset-DFS + the fresh slot-availability CNF solved by kissat
  4.0.4): **agreement on all 1,252, zero disagreements**.
- Measured cost **110.7 CPU-seconds** of the pre-registered **6.0 CPU-hour**
  cap (single thread, `nice -n 15`, BLAS/OMP pinned to 1,
  `PYTHONDONTWRITEBYTECODE=1`). Nothing was budget-bound; the negative is
  COMPLETE, not truncated.
- Checkpointing: one fsynced line per decided instance
  (`scan_U_checkpoint.jsonl` 428, `scan_V_checkpoint.jsonl` 426,
  `scan_Wfac_checkpoint.jsonl` 398).

## Controls (both directions, all in-run, all PASSED)

ACCEPT — A1: 729/729 Brent identities over Z in `fmpz` AND pure int, all
factors ternary; A1b frozen landscape row reproduced through the same code
path (d = 14/14/13, floors F/F/F, DFS states 1152/1344/460, total_lb 58);
A2 printed-SLP recount by exact expansion = 16/16/28 = 60 with every stage
reproducing its target rows (23 U rows, 23 V rows, all 9 C-entries exactly);
A3 extension-positive control on an INDEPENDENT target (tau(sun56-V) + Sun's
two aux classes at T = 13) positive in **three** instruments — DFS, kissat,
CaDiCaL 3.0.1; A4 transposition preconditions audited (23/23 nonzero rows on
all three sides, rank(Wfac) = 9 exactly over Q).

REJECT — R1 one-addition-deleted plant on the verified 16-gate printed-left
carrier: verifier REJECTS (one tau class uncovered); R2 operand-sign
perturbation of a consumed gate: verifier REJECTS by exact recomputation;
R3 weaker-T negative (tau(Wfac) at T = 12) INFEASIBLE with a fresh
DRAT→LRAT certificate accepted by both pinned checkers; R4 planted wrong
d-count assertion fires (rc 1).

## Certificates

**37 DRAT→LRAT certificates**, every one accepted by BOTH independently
pinned checkers (`drat-trim` `111b0405…` `s VERIFIED`, `lrat-check`
`b4bdebfc…` rc 0) — the three floor instances (U@14, V@14, Wfac@13), the R3
weaker-T instance, and for each scan the lexicographically first infeasible
instance plus a seed-20260902 sample of 10. All 37 re-verified on replay by
the independent audit after the run.

## Independent post-run audit — 7/7 blocks PASS (`postcompute_audit.json`)

Third-path universe re-enumeration with a different loop shape (all 6 hashes
match the pins), checkpoint completeness (428/426/398 with 0/0/8 admissions),
24 boundary+random DFS re-runs matching their checkpointed verdicts, the
14-gate Wfac witness re-verified by exact expansion, the printed 16/16/28 = 60
recount re-derived independently, all 37 certificates replayed through both
checkers, and the verdict arithmetic plus Brent anchors re-checked.

## Frontier statement (bounded, precise)

For the `stapleton60` factor blocks in their fixed orientation (sigma^0,
all-monomial data triple), in the three-stage linear SLP model (inputs free;
gate = x±y; sign changes and copies free; output counted with the
transposition bound): **the minimum total is exactly 60**, attained by
Stapleton's own published circuit. Two consequences worth stating plainly:

1. The frozen certified LB **58 was not tight** — both side floors' d+1
   bounds (15 and 15) are unattainable; the true side costs are 16 and 16.
   Nothing frozen is contradicted: 58 was and remains a valid lower bound.
2. This is a *negative* frontier move on the upper-bound question (no 58 or
   59 circuit exists here) and simultaneously a *positive* optimality
   certificate for a published scheme — the same shape as the frozen sun56
   result, now at the other end of the ladder.

**Rule-7 scope.** Decided exactly: "C(U) <= 15?" over the pinned 428-element
aux-1 universe at T = 15, "C(V) <= 15?" over the pinned 426-element universe
at T = 15, and "C(Wfac) <= 14?" over the pinned 398-element universe at
T = 14, for this one orientation and data triple, under this model. NOT
searched: other stapleton60 data triples or sigma classes (the frozen census
gives their lower bounds only), other decompositions, non-ternary alphabets,
GL(3,Q)/GL(3,Z) sandwiches, other addition-count conventions, and the open
`mws59` 58-vs-59 residue. `paper55`'s 55 remains the record, untouched.

## Instrument defects

None new in this run. Every control and certificate passed on the first
execution; no instrument was rebuilt mid-campaign, no verdict-relevant code
was changed after init. Carried forward from session 10 and honoured here:
the frozen session-3 `transpose_check.py` is void, so **this campaign used no
transposition-derived witness** — the C(Wfac) = 14 upper bound is a
synthesized explicit 14-gate circuit verified by exact expansion, and the
output-stage 28 is the printed stage re-expanded, not a transposed artifact.

## exactly-one verdict

**FROZEN-NEGATIVE** — no 58- or 59-addition circuit exists for
`stapleton60`'s fixed orientation; the searched space is the three pinned
auxiliary universes (428 + 426 + 398 = 1,252 dual-instrument instances, zero
disagreement, 110.7 CPU-s of a 6.0 CPU-h budget); consequently the
fixed-orientation minimum is **exactly 60** (LB = UB = 60) and the published
scheme is certified optimal for its own orientation.
