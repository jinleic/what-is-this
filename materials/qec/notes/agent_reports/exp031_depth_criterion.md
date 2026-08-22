NEGATIVE - the depth criterion and the two-value law are FALSIFIED in the translation-invariant class by 3 weight-9 PBB members; they hold on all 208 decided instances. [Rewritten by parent from the regenerated artifact per FR-014; the original POSITIVE report filed certified refutations as 'undecided'.]

# EXP-031: one-ancilla depth criterion - corrected report (2026-08-12)

## What was tested
Criterion: minimum one-ancilla syndrome-extraction depth = w_max + [parity
obstruction], where the obstruction is exact infeasibility of the depth-w_max
incidence-colouring + Lemma-C1 CSP (independent PySAT cadical195 implementation);
ground truth is cpsat_schedule with INFEASIBLE certificates at every smaller
depth. Schedule classes are explicit: BB/PBB tiers translation-invariant (TI),
random tier unrestricted.

## Headline result (corrected)
- **Verdict: NEGATIVE.** criterion FALSIFIED: 3 instance(s) violate the proposed depth criterion (certified in-class refutations count even when the true minimum remains undecided)
- The three counterexamples are pbb144-01, pbb144-03, pbb144-11 (all w_max=9):
  CP-SAT certifies INFEASIBLE at EVERY depth 9,10,11,12,13 in the TI class.
  The criterion predicted depth 10; a certified infeasibility at the predicted
  depth refutes the prediction regardless of the true minimum (FR-014).
  Consequently these members admit NO translation-invariant schedule at any
  depth <= 13 - the two-value law {w_max, w_max+1} is REFUTED in the TI class.
- On every DECIDED instance the criterion is exact: 208/208 matches,
  0 mismatches among decided; reduction law 0 mismatches; two-value law holds
  on all decided instances.
  Tiers: 47 BB (TI), 11/14 PBB (TI), 150 random commuting sets (unrestricted).
- Scope-corrected statement that the data supports: *within the tested tiers,
  the criterion predicts the certified minimum wherever a minimum exists in
  class; for 3 of 14 PBB members the TI class itself is empty through depth 13,
  so any TI-scoped depth statement must exclude them.*
- Whether the two-value law survives UNRESTRICTED for those 3 members is
  EXP-033 (in flight at report time; see its artifact for the outcome).

## Weak closed-form predicate ('exists anticommuting overlap')
TP 51 / FP 157 / TN 3 / FN 0
(false-positive rate 0.98125); exact on the 47-code BB family only.
It is NOT the criterion; the exact CSP is.

## Prior art
ASC (arXiv:2603.21499, 2026-03-23) certifies no-depth-6 for IBM's BB codes via
SMT but states only the Tanner max-degree lower bound analytically - it does
not explain the 6->7 gap nor treat mixed checks. Our contribution is the
analytic parity criterion + the class-scoped falsification boundary above.

## Bookkeeping
- Attempted 211; decided 208; undecided 3
  (all three undecided are the certified TI refutations, n_refuted_undecided=3).
- Canonical artifact: results/processed/exp031_depth_criterion.json
  sha256 3c52b17a7cac98bed5813d93586fa1397028469e886ec05c666088fabc89af6e
- wall_s 76.231
- Regression tests: tests/test_exp031_verdict.py (8 pass) - FR-013 polarity +
  FR-014 refutation classification.
- Reproduce: PYTHONPATH=src .venv/bin/python experiments/exp031_depth_criterion.py
