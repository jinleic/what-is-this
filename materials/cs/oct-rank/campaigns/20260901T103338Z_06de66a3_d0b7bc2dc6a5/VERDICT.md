# N1 VERDICT.md — exact CP-completion campaign (OctRankNext)

Campaign: `campaigns/20260901T103338Z_06de66a3_d0b7bc2dc6a5/` — gate
`n1_exact_cp_completion`. Prereg committed BEFORE any compute:
commit `d6c7e44` (prereg sha256 `d569057e…f82`, byte-identical copy as
`pre_statement.md` in this dir, PROVENANCE.md records the copy discipline).

## Headline

**FAILURE TO CERTIFY — the N1 exact-CP-completion attack produced no
admitted candidate and no rank-13 witness. rank((L_1, L_i, L_j)) stays
OPEN in {13, 14}. By the pre-registered vocabulary: "Absence of a 13-witness
is NOT evidence for 14; absence of an impossibility argument is NOT evidence
for 13." The published window 18 <= R_R(T_O) <= 25 is untouched; no
published work is refuted.**

## What was fixed before compute (prereg, decision-complete)

- Rank fixed: r = 13 only. Field/ring: exact fmpq on every load-bearing
  number; float64 only for polish, slice choice, diagnostics; final
  comparisons outward-arb re-verified (frozen discipline).
- Normalization: geometric-mean column balancing between rounds (frozen
  gauge move, image preserved exactly). No symmetry ansatz imposed; the
  blockdiag structure is a recorded anchor, not a constraint.
- Search cap: 40 seeds (32 random base-20260901, 7 tau-merge, 1
  frozen-restart), 3 TRF rounds x 5000 nfev, admission ||g0||_max <= 1e-7
  (exact), exact Newton (3 iters max), frozen 15-rung ladder; 2 h CPU cap.
- Kill criteria: E4 fires on zero admissions or non-opening window; verdict
  mapping fixed in prereg section 7/9.
- Instrument: byte-copy of the frozen Route-F certification core
  (n1_instrument_base.py sha256 781e51ec1a72cfb8… == frozen rf_krawczyk.py).

## Controls (run BEFORE the sweep; behavior contract, prereg section 4)

| control | role | outcome (frozen n1_results.json) |
|---|---|---|
| C-TAU7 | accept known-valid (tau, r=7, frozen certificate) | CONTAINMENT at rho=1e-3 — instrument certifies a true case |
| C-TAU6 | reject globally infeasible (tau, r=6; true rank 7 machine-verified) | NO-CONTAINMENT-ANY-RUNG — no false certificate |
| P-POS | accept synthetic known-true (rank-13 plant, S-only 1e-8) | CONTAINMENT at rho=1e-4 |
| P-WRNG | reject ONE-COORDINATE perturbation of the target ([0,0,0]=2; the assignment's named near-miss plant) | NO-CONTAINMENT at every rung; complete certified exclusion rho<=1e-4 |
| behavior_ok | quarantine criterion | true (never triggered) |

The instrument accepts a known-valid decomposition (tau-7, P-pos) and
rejects a one-coordinate perturbation (P-wrong, tau-6) — the assignment's
control contract, MACHINE-VERIFIED through the frozen code path.

## Main sweep (all 40 seeds completed; none skipped; 3769.7 s CPU < 2 h cap)

| class | n | best rel (final) | exact ||g0||_max (range) | admitted |
|---|---|---|---|---|
| R (random) | 32 | 1.414e-06 | ~1e-5 .. 5e-4 | 0 |
| M (tau-merge) | 7 | 2.472e-06 | ~1e-5 .. 3e-4 | 0 |
| F (frozen-restart) | 1 | 6.886e-06 | ~4.2e-5 | 0 |

- Best seed: R18 (rel 1.414e-06). The full 40-row table with exact gmax
  fractions is in n1_results.json (frozen).
- Class pattern matches the frozen Route-F stall: TRF bottoms at
  rel ~1e-6..4e-5 for EVERY pip — including the 7 merge seeds that
  initialize at the certified 14-term structure's neighbourhood and the
  frozen-candidate restart probing whether a deeper root exists in the same
  basin (the F seed's final rel 6.9e-06 is WORSE than the frozen candidate's
  8.8e-06 polished rel was — no deeper basin was reached; note the F probe
  cannot contradict the frozen certified exclusion, which covered only
  rho in [1e-6, 1e-12] boxes).
- With 0 admissions: exact Newton never engaged; the ladder never ran on a
  main candidate; E4 kill fired exactly as pre-registered.

## Adjudication (prereg section 7 vocabulary)

Verdict: **FAILURE TO CERTIFY**. No rank claim in either direction.
Both frozen absence sentences apply verbatim. N1 does NOT make N2 logically
unnecessary: no witness was produced, so the second preregistration clause
(implication fires only on a witness) does NOT trigger — N2 (upward
Gröbner/elimination, infeasibility direction) proceeds as a distinct
campaign under its own committed prereg, per the assignment sequencing.

## Defects disclosed (rule 5)

1. Run 1 crashed pre-controls (tauT target builder emitted 4-vectors instead
   of scalars — a code defect, not a convention error). Fixed; the broken
   run's output is preserved as n1_run.out.run1_broken; ZERO certification
   or control data was produced by run 1 (crash inside the first certify()
   call). Run 2 (of record) re-ran everything from scratch.
2. No other defects; no post-hoc parameter changes (kill fired on
   pre-registered criteria; admission threshold never adjusted).

## Cost (in-process accounting, rule 17e)

time.process_time() total 3769.7 s CPU, single-threaded, nice -n 10.
No ps-derived numbers anywhere. All scripts exit cleanly; no daemonized
processes.

## Rule-7 scope sentence

Covered: the fixed conjugated Route-F tensor (blockdiag(tau,tau), proven
rank-equivalent to (L_1,L_i,L_j) by two frozen entrywise diagonal
identities) at rank 13, from exactly the 40 declared seeds through the
declared 3-round TRF polish, exact fmpq residual admission at 1e-7, the
declared exact-Newton protocol (never engaged — no admission), the frozen
15-rung ladder (ran only inside the four controls), and the four declared
controls through the byte-frozen instrument. NOT covered, hence not excluded
by anything here: any rank-13 witness anywhere (no witness was produced);
any rank >= 14 statement; seeds/init schemes beyond the 40; admission
thresholds other than 1e-7; anything off the declared polish budget; the
commuting-extension model; border rank; complex decompositions; other
octonion triples; and anything touching R_R(T_O) itself
(18 <= R_R(T_O) <= 25 stays untouched; no published work is refuted).
