# Independent UC analytic referee report — 2026-08-26

Invocation: ephemeral `codex exec`, model `gpt-5.5`, high reasoning,
read-only sandbox.  Scope was restricted to `paper/main.tex`, `AUDIT.md`,
`PROOF.md`, `bridge_uc.py`, `thmB3_proof.py`, `margin_lemma.py`, and
`lemma_rh_proof.py`.  The reviewer was instructed not to trust the existing
audit, to try to falsify the human proof, and not to inspect long trace
arithmetic.

## Verdict on the pre-repair text

- Blocking findings: none.
- Mathematical chain: no counterexample or reversed implication found.
- Two major exposition findings and two minor precision findings, all repaired
  in the current sources.

## Findings and resolutions

1. `PROOF.md` multiplied the diagonal defect by entropy twice:
   `||psi||^2 = rh*h`, although that document defined `rh` as the complete
   diagonal defect.  It now states the correct identity
   `||psi||^2 = W(x,x)/ln(2) = rh = rho*h`.
2. The manuscript stated equality between the coupling and pair-orbit
   optimizations after proving only the fold direction.  It now explicitly
   unfolds every pair-orbit law: off-diagonal mass is split equally between
   `(p,q)` and `(q,p)`, diagonal mass remains diagonal, both marginals are
   `mu_nu`, and cost/folding are preserved.
3. The scale-free paragraph in `margin_lemma.py` stated equality with
   `kappa_t` for arbitrary sink collapse.  It now proves the correct general
   statement.  If `epsilon` is non-sink mass, then
   `sigma^2 <= epsilon*rho_max*L`; for limiting mean `m<=t`,
   `liminf Lambda >= 2(1-alpha)(1-m)-1 >= kappa_t`.  Equality is reserved for
   the displayed mean-pinned zero-cost path.  The manuscript carries the same
   correction.
4. The direct extreme-point perturbation argument now says explicitly that the
   three chosen points lie in the topological support, so their disjoint Borel
   neighborhoods have positive measure.  Both the manuscript and
   `thmB3_proof.py` were corrected.

## Independently checked transitions

The referee checked the exact `s*` equivalence; symmetrization, folding, and
reverse unfolding; compact optimal-coupling attainment; equality of the
infinite-dimensional and pair-orbit optimizations; compactness and continuity
on `P(Delta)`; fixed-mean concavity; Bauer and one-moment extreme-point steps;
`Phi_exact-Phi_rel=(1-alpha)(S^2-E)>=0`; the Bochner Margin Lemma and `L=0`;
the `rh` inequality and endpoint signs; five-parameter domain and orbit-swap
coverage; the sequential Bernoulli entropy argument; and the strict/non-strict
boundary at the certified target.

Post-repair executable controls:

- `python -I -B uc/bridge_uc.py`: `BRIDGE CHECK PASS`;
- `python -I -B uc/thmB3_proof.py`: all exact/cited setup checks and three
  sampled controls pass;
- `python -I -B uc/margin_lemma.py`: `ALL PASS` (its final historical sentence
  remains a pre-Campaign-I chronology marker, not current status);
- the 35-page manuscript rebuilds without undefined citations or references.

## Residual assumptions

This is independent-model review evidence, not journal peer review or formal
proof.  The named standard functional-analysis theorems, the separately
reviewed decomposition/reduction identities, Arb primitives, and the finite
certificate remain outside this review's human-only scope.
