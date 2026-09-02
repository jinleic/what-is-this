# Campaign report — gate `paper55-sigma12-total55` (run 20260902T015158Z_c0bea75d_761cb990fc40)

## Question (prereg sha `10b3f32f…`, commit `bdf3a13`)

Decide the minimum of the `paper55` **sigma^1** and **sigma^2** orientations,
each frozen at **[55, 57]** (certified LB 55; best known circuits 57).

## Answer — verdict FROZEN-CERTIFIED: both orientations are EXACTLY 55

| orientation | left | right | output | total | LB/UB after |
|---|---|---|---|---|---|
| sigma^1 = (V, T(W), T(U)) | 14 | 14 | **27** | **55** | 55 / 55 |
| sigma^2 = (T(W), U, T(V)) | 14 | 13 | **28** | **55** | 55 / 55 |

Both "best known 57" rows are superseded by 2. **This is the record VALUE at two
further orientations of the same tensor, not a new record** — `paper55`
sigma^0's 55 remains the record, and nothing about sub-55 is claimed.

## Why no search was needed for the side costs (the reduction)

The frozen gate-B results make all three sigma^0 side costs exact:
C(U) = 13, C(V) = 14, C(Wfac) = 14. The campaign's key structural step is the
**proved free-relabeling symmetry**: applying a fixed permutation of the input
coordinates carries any circuit to a circuit of the same cost (rename the free
input wires), so `C(T(X)) = C(X)` for the per-product transpose `T`. Re-checked
numerically in-run (control A4: d and the floor verdict of `T(X)` equal those of
`X` for X ∈ {U, V, W}). Hence:

- sigma^1: left = V (14 exact), right = T(W) (14 exact), output-factor = T(U)
  (13 exact) ⇒ output >= 27 ⇒ total >= 55 — the frozen LB, now with every side
  cost exact rather than floor+1;
- sigma^2: left = T(W) (14), right = U (13), output-factor = T(V) (14) ⇒
  output >= 28 ⇒ total >= 55.

So the only open quantity per orientation was whether its output stage attains
its bound. The campaign closed both by constructing and verifying explicit
circuits.

## What was constructed, and how each piece is verified

1. **Wfac circuit synthesized (not transposition-derived).** The aux-1 dual
   instrument scanned the hash-pinned universe |AU(Wfac)| = **391**
   (sha256 `1a3e66e9…`, two independent enumeration loops agreeing) at
   T = d+1 = 14; instance 307 (aux `(0,0,-1,1,0,-1,-1,0,1)`) admitted, its
   schedule was replayed as an explicit **14-gate circuit** and verified by
   exact expansion against all 23 Wfac rows (`witness_Wfac_14gates.json`).
   Per the session-10 D1 lesson, no transposition-derived Wfac witness is used
   anywhere.
2. **Left/right stages** are the paper's printed circuits under free input
   renaming, each re-verified in-run by exact expansion: `RIGHT_SLP` 14 gates →
   23 V rows; `LEFT_SLP` 13 gates → 23 U rows; the synthesized 14-gate Wfac
   circuit under the T relabeling → 23 T(W) rows.
3. **Output stages** were built by **reverse-mode transposition** implemented
   here as explicit adjoint accumulations, with additions counted exactly (each
   accumulation into a non-empty adjoint costs 1): **27** additions for
   sigma^1 (from the 13-gate U circuit) and **28** for sigma^2 (from the
   14-gate V circuit). The construction is not evidence: each constructed stage
   was checked to reproduce its orientation's W-combination map entry by entry
   (23 × 9 exact equalities per orientation).
4. **End-to-end.** Each assembled scheme was expanded symbolically over all
   **81 monomials A_i·B_j** and compared against the bilinear map defined by
   that orientation's frozen factor triple: **9/9 outputs exact, zero
   mismatches, both orientations** (`assembled_sigma1.json`,
   `assembled_sigma2.json`). All three orientations' triples were also
   re-verified 729/729 Brent over Z in `fmpz` AND pure int.

## Controls (both directions, all in-run, all passed)

ACCEPT — A1: 729/729 Brent for sigma^0, sigma^1, sigma^2 in two arithmetics;
A2: the printed SLP re-verified and recounted (13 + 14 + 28 = **55**, each stage
reproducing its targets, output stage reproducing the W columns exactly);
A3: frozen landscape rows reproduced through the same code path — d/state
counts (12,33), (13,116), (13,66) for sigma^0 and their cyclic images for
sigma^1/sigma^2, floors impossible everywhere; A4: the relabeling symmetry
re-checked numerically; A5: transposition preconditions audited per orientation
(23/23 nonzero rows, rank of the output factor = 9 exactly over Q); plus a
tri-instrument SAT control (DFS + kissat 4.0.4 + CaDiCaL 3.0.1) on an
independent target.

REJECT — R1 one-addition-deleted plant on the printed left circuit: REJECTED
(uncovered target class); R2 operand-sign perturbation: REJECTED; R4 the frozen
Wfac floor at T = 13 re-certified INFEASIBLE in both instruments with a fresh
DRAT→LRAT certificate accepted by both pinned checkers (rc 0/0); R5 planted
wrong d-count assertion fired.

## Independent post-run audit — 6/6 blocks PASS (`postcompute_audit.json`)

Floors re-derived; the relabeling symmetry re-checked; the synthesized Wfac
witness re-verified; **both assemblies re-derived from scratch and re-verified
end to end (55 each, 9/9 outputs)**; the certificate replayed through both
checkers; scan completeness and verdict consistency.

## Instrument defects (disclosed; both caught by pre-registered checks, no verdict was affected)

- **D1 (attempt 1, fixed):** the first assembly attempt relabeled circuits by
  rewriting their stored values, which the exact-expansion verifier then
  recomputed from unpermuted inputs — the verifier correctly REFUSED the stage
  (`sigma^1: stage verification failed R=False`). Fixed by threading the input
  permutation into the verifier/evaluator instead of rewriting values.
- **D2 (attempt 2, fixed):** with a permuted input assignment, the transposed
  stage's 9 outputs are indexed by input WIRES, so they emerged in permuted
  coordinate order; the pre-registered output-map equality check caught it
  (`OUTPUT-MAP-MISMATCH` with the counts already correct at 27/28 and totals
  55/55). Fixed by reporting the transposed outputs in coordinate order —
  re-indexing output wires is a free relabeling. The attempt-2 INCONCLUSIVE
  verdict is preserved verbatim (`verdict_attempt1_inconclusive.json`,
  `runner_output_attempt1.json`); no claim rests on either attempt.
- The Wfac scan ran in each attempt (3 complete passes of 391 instances, 1,173
  checkpoint rows, 3 admissions — the same instance 307 each time), all
  dual-instrument agreeing; only the final pass backs the verdict.

## Cost

Total measured **≈ 81 CPU-seconds** for the verdict-bearing pass (391
dual-instrument instances plus arithmetic) of the pre-registered **2.0 CPU-hour**
cap; no budget was extended, nothing was sampled.

## Frontier statement (bounded, precise)

In the three-stage linear SLP model (inputs free; gate = x±y; sign changes and
copies free; output counted with the transposition bound), for the `paper55`
factor blocks: **the sigma^1 and sigma^2 orientations each have minimum exactly
55 additions**, with explicit verified circuits (14/14/27 and 14/13/28). Their
frozen intervals [55, 57] collapse to {55}. No claim is made about other
decompositions, other data triples or sandwiches, non-ternary alphabets,
GL(3,ℚ)/GL(3,ℤ), other counting conventions, or the global sub-55 question.

## Exact fixed-orientation ladder after sessions 9–13

| decomposition / orientation | published | certified exact minimum |
|---|---|---|
| `paper55` sigma^0 | 55 | **55** (record) |
| `paper55` sigma^1 | — (57 best known) | **55** |
| `paper55` sigma^2 | — (57 best known) | **55** |
| `sun56` sigma^0 | 56 | **56** (published optimal) |
| `mws59` sigma^0 | 59 | **58** (published +1) |
| `stapleton60` sigma^0 | 60 | **60** (published optimal) |

## exactly-one verdict

**FROZEN-CERTIFIED** — both `paper55` sigma^1 and sigma^2 attain **55**, each
with an explicit circuit verified by exact expansion stage-by-stage and end to
end over all 81 monomials; LB = UB = 55 for both orientations; the record value
now holds at three orientations of the same tensor, and the record itself is
unchanged.
