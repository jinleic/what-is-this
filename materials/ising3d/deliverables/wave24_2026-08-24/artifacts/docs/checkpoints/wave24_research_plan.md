# Wave 24 research plan

**Status:** frozen before implementation in the Wave-24 session; materialized here during integration.
**Date:** 2026-08-24.

## Goal

Push three independent exact fronts rather than add another same-point computation:

1. turn the wave-22 qualitative finite full-spectrum exceptional sets into a graph-uniform
   algebraic size theorem;
2. extract a genuinely positive structural theorem from exact checkerboard decimation while
   locating the positivity/topology obstruction honestly;
3. test a non-termwise Callen subsystem by combining selected rows before any support closure.

A result lands only with a producer, a clean-room verifier, a proof with explicit scope, current
source hashes, a machine artifact, adversarial review, ledger/checkpoint integration, and a
validated evidence bundle. Exact integer/rational/algebraic arithmetic only. `K_c=0.221654626`
remains comparison-only.

## Research lanes surveyed first

- **Spectral/free-fermion:** matchgate identities, pure sub-Pfaffian tensors, hidden free-fermion
  circuits, and algebraic projection/degree theorems.
- **Correlation bootstrap:** lattice spin-flip/Callen equations, symmetry-reduced row systems,
  Fourier blocks, Krawtchouk transforms, and convergence results that must not be mistaken for
  finite closure.
- **Surface/decimation:** normal-factor-graph duality, positive even-subgraph expansions,
  higher-spin interactions, and topological-sector caveats.
- **W-law algebra:** inclusion matrices, Johnson/Bier transforms, Koszul/incidence complexes,
  and the exact 17/133 residual cores.

Primary inputs selected for direct use were the 2025 projective degree theorem
[arXiv:2511.07639v1](https://arxiv.org/abs/2511.07639v1), the NIST DLMF cosine-hyperbolic product
[DLMF 4.36.E2](https://dlmf.nist.gov/4.36.E2), and the lattice-bootstrap papers
[arXiv:2206.12538](https://arxiv.org/abs/2206.12538) and
[arXiv:2309.01016](https://arxiv.org/abs/2309.01016). Literature suggested techniques; no external
numerical value was used as theorem evidence.

## Ranked candidates and decisions

### A. Exceptional-scheme degree — execute

- Polynomialize a scalar-equivalent isotropic transfer matrix.
- Bound every coefficient-incidence equation in the `n+2` variables `(t,c,u_1,...,u_n)`.
- Combine the existing complex finite projection with an ambient-dimension degree theorem.
- **Success:** an explicit bound for every connected simple nonpath graph, with paths retained as
  a cofinite control.
- **Failure signature:** saturation/projective components or path scope invalidate the finite image.

### B. Even-star/matchgate structure — execute, high risk

- Expand the eliminated star before taking its logarithm.
- Test whether the positive even Walsh tensor is a principal-sub-Pfaffian matchgate.
- Separate ordinary flattening rank, nonnegative all-leg CP rank, and the trivial operation of
  restoring the eliminated centre.
- Search for an all-degree sign law rather than stop at `c2,c4,c6`.
- **Success:** self-contained analytic theorem plus exact finite Pfaffian/rank/cycle controls.
- **Failure signature:** a violated matchgate identity, sign change, or a nonnegative rectangle
  covering more than one even state.

### C. Selected-row Callen Fourier block — execute, high risk

- Use a far marked spin so pivot rows have no contact terms.
- Form the proper-even by odd neighbour-cube leakage matrix and diagonalize it by Walsh characters.
- Fold the cubic case by inversion and construct explicit row combinations.
- **Success:** a finite exact compressed relation or a theorem-grade obstruction naming exactly why
  the pair row cannot close.
- **Failure signature:** no leakage kernel, wrong rank under folding, or support/contact growth with
  volume.

### D. W-law block core — investigate, defer implementation

The exact cores decode as complement-packing/empty-rung boundary strata, but no uniform differential,
block determinant, or contracting homotopy was proved. Another leaf order or fitted determinant was
forbidden. The next experiment is a support-defect filtration, not an `L=10` closure.

## Frozen file ownership

- `e244`: octahedral matchgate/sign producer, verifier, artifact, proof.
- `e245`: exceptional-degree producer, verifier, artifact, proof.
- `e246`: Callen Walsh producer, verifier, artifact, proof.
- Shared summaries, append-only ledger, checkpoints, reports, and bundle: lead-only integration.
- The pre-existing `e233_trace_resultant_positivity.py` draft stays unlanded and excluded.

## Verification contract

- Producers fail before writing JSON if any mathematical, provenance, scope, CPU, or RSS gate fails.
- Verifiers import neither their producer nor implementation helpers that would duplicate the same
  failure mode.
- Universal finite tables have complete row coverage; source hashes are hard predicates.
- Review findings cause source/artifact regeneration and a fresh verifier run.
- The full suite is not required under severe shared-host load; every changed theorem surface must
  pass its targeted producer/verifier and independent review.

## Outcome against plan

All three selected fronts landed. Review and completion audit corrected six load-bearing issues
before finalization: retained-spin normalization, path/cofinite scope, inherited dependency
provenance, the raw/folded Callen deleted-row rank interpretation, e245 RSS-method provenance,
and registration of every Wave-24 external source. The W-law lane produced a sharper next
hypothesis but no claim. No thermodynamic observable or critical endpoint moved.
