# Wave 25 research plan

**Status:** frozen before implementation.
**Date:** 2026-08-24.

## Goal

Push three independent exact fronts under a strict single-process resource wall:

1. expose the first global obstruction created when locally Gaussian Ising-star tensors are
   contracted through nonplanar three-dimensional geometry;
2. attack the open `2x4` full-spectrum exceptional scheme with the first viable trace-nine
   replica-column engine, correcting the tempting but underdetermined trace-eight plan;
3. turn the observed `17=16+1` and `133=128+5` W-law cores into a proved support-defect
   decomposition and seek a block/Schur determinant mechanism rather than another leaf order.

A theorem lands only with an exact producer, a clean-room verifier, a scoped proof, a machine
artifact, current source hashes, adversarial review, ledger/checkpoint integration, and a validated
evidence bundle. A failed high-risk route lands only as a reproducible negative result with the
first exact failing gate. No thermodynamic solution or critical-point movement is inferred from a
finite tensor, finite layer, or resource-bounded computation.

## Research surveyed before the freeze

- **Recent fermionic tensor networks.** Arbitrary-geometry fermionic contraction makes ordering and
  swap structure explicit ([arXiv:2410.02215v3](https://arxiv.org/abs/2410.02215v3)); a minimal
  tensor beyond free fermions isolates the interacting deformation
  ([arXiv:2412.04216v2](https://arxiv.org/abs/2412.04216v2)); parity-preserving circuits admit a
  Gaussian-plus-universal-quartic decomposition
  ([arXiv:2504.19317v1](https://arxiv.org/abs/2504.19317v1)); and Quon simulation frames matchgate
  complexity topologically ([arXiv:2505.07804v2](https://arxiv.org/abs/2505.07804v2)). These papers
  motivate the crossing-defect experiment; no external numerical claim is used as evidence.
- **Recent 3D Ising tensor work.** The 2026 interface study reduces a highly asymmetric 3D problem
  to an effective 2D tensor network rather than claiming an exact bulk solution
  ([arXiv:2601.07829v1](https://arxiv.org/abs/2601.07829v1)). This reinforces the need to preserve
  finite-geometry scope.
- **Moment/localizer route.** Shifted Stieltjes matrices are the correct semialgebraic test for mode
  roots `u_i>=4`, but inequalities do not lower the complex incidence dimension. For eight modes,
  `Q(0),Q(2),Q(3)` leave five coefficient variables; `F4` through `F8` only form the coefficient
  ideal over `Q(q)`, and `F9` is required to cut the free coupling. Any trace-eight finite-root
  claim is forbidden.
- **W-law route.** The stored exact cores are boundary packings in the rung-defect filtration, but
  equality between the maximal leaf residual and that boundary is known only at `L=8,9`. Raw
  singleton peeling, the `#a<=2` truncation, and an unsigned naive Koszul differential are closed.

## Executed lanes and stop rules

### A. Boundary matchgate topology — execute first

- Checkerboard-eliminate one color of the smallest nonplanar open simple-cubic box, `2x3x3`.
- Retain a face boundary, average all other retained-color spins, and compute its exact Walsh
  signature as a polynomial in `v=tanh K`.
- Reconstruct the candidate skew pair matrix and test every boundary ordering against all
  principal-Pfaffian identities. Include a planar alternating-cycle control.
- **Breakthrough target:** one explicit identity residual with a fixed nonzero sign for every
  `0<v<1`, proving a genuine global crossing defect despite local matchgate tensors.
- **Scope boundary:** failure of one boundary tensor to be one matchgate is not a Pfaffian-number
  lower bound and does not exclude spin-structure sums or hidden auxiliaries.

### B. Open `2x4` trace-nine exceptional scheme — execute aggressively, resource-gated

- Prove and implement the replica-column identity for `tr(B(q)^k)` on open `2xL` without dense
  symbolic matrix powers. Validate it bit-for-bit against the landed `2x3` traces.
- For `2x4`, generate exact `H1,...,H9`; the largest local state space is `4^9=262144` and degree is
  at most `162`. Use compact modular/CRT or packed coefficients; never two Python-integer arrays of
  all `42.7` million slots.
- Build the degree-eight mode polynomial. Use `F4,...,F8` for the five-variable coefficient ideal;
  only if its leading-form/quotient and memory preflights pass, form the `F9` multiplication norm.
- Dispose candidate physical branches with Hermite signature plus the shifted localizer
  `H1=[m_{i+j+1}]` for `y_i=u_i-4`.
- **Breakthrough target:** prove the open `2x4` physical exceptional set empty.
- **Stop rule:** one nice/low-priority process, process CPU at most `900 s` per stage, task RSS below
  `2 GiB`; abort before a predicted `>1 GiB` materialization, nonregular quotient, or unverified
  dimension jump. Record the exact first blocked gate rather than narrowing the theorem silently.

### C. W-law support-defect block — execute in parallel conceptually, serially on CPU

- Decode every existing `17` and `133` core coordinate into `(T,B,E,D)` rung data and prove
  `E-D=L-2m` and `E+D>=|L-2m|`, with equality exactly on the packing boundary.
- Certify the finite decompositions `17=16+1` and `133=128+5` by Burnside counts and direct artifact
  coverage, including particle-hole duality of `(9,4)` and `(9,5)`.
- Partition the literal core matrices by packing/boundary stratum and seek exact row/column
  operations or a Schur complement explaining `-2^32*3` and `-2^288*3^5`.
- **Breakthrough target:** a uniform block lemma that implies full column rank beyond the two stored
  sizes.
- **Fallback result:** a proved finite support-defect decomposition plus an exact obstruction to the
  proposed Schur recurrence, with no all-`L` promotion.
- **Stop rule:** use only stored `L<=9` annihilators; do not launch an `L=10` closure or fit a
  recurrence from two determinants.

### D. Contact-inclusive Callen block — reserve lane

If one executed lane fails early, form the cubic-symmetry quotient of all degree-six contact rows on
`C3^3` and decide exactly whether any combination leaves only normalization and nearest-neighbour
pair columns. This graph is intentionally non-bipartite; it is valid for Callen identities but never
for checkerboard elimination. A result is finite and contact-inclusive only.

## Resource and verification contract

- Live pre-freeze load was `10.95 / 10.88 / 6.08` on 28 cores with ten users. Heavy work remains
  single-process and low priority; no full-suite launch.
- Producers use process CPU and Darwin `mach_task_basic_info.resident_size_max`, fail before writing
  on a resource or semantic gate, and record the measurement method.
- Verifiers do not import producers. Universal tables require complete coverage and source hashes are
  pass predicates.
- Reviews must attack dimension counts, normalization, ordering/sign conventions, positive controls,
  and scope. Review findings trigger regeneration and a fresh targeted run.
- The pre-existing unlanded `e233_trace_resultant_positivity.py` remains excluded. Wave 25 begins at
  `e247` and ledger row `H658`.

## Outcome against plan

All three primary lanes produced exact results, and the reserve Callen lane was activated after the
spectral projection remained open.

- **Lane A landed stronger than its rational-point fallback.** Both checkerboard colors on open
  `2x3x3` give the same boundary signature, and four exact factor classes prove failure of every
  ordinary unsigned leg order throughout `0<v<1`. The planar control passes.
- **Lane B moved the bottleneck but did not empty the scheme.** Packed trace generation through
  `k=9` succeeded in 262,144 states. Independent Lucas and physical-mode reconstruction plus a
  filtered/projective Bezout lemma prove the `F4`-`F8` coefficient scheme has length 96.
  The explicit lifted basis, `96x96` `F9` norm, and physical branch disposition remain open.
  An optional actual-quotient probe reached its outer 900-second background wall and wrote no
  artifact.
- **Lane C landed the finite mechanism and falsified its simplest recurrence.** The 17/133 cores
  are exact packing-boundary families; L9 is four 29-shells plus the signed inherited L8 core.
  Smith factor 6 kills the proposed eight-copy/dyadic recurrence. Nothing is promoted to L>=10.
- **Reserve lane D landed an all-positive-coupling no-go.** All 864 radius-one contact rows on
  `C3^3`, including the size-five neighbour-contact reduction of e246's `U=N` row, have a positive
  nonallowed cubic-orbit minor on `0<v<1`.

All four producers and all four clean-room verifiers passed. Four adversarial reviews returned final
PASS after repairs to ordinary-order scope, actual trace-degree and overflow predicates, independent
Lucas/mode linkage, full semantic hashes, dependency locks/runtime versions, and Callen
pivot/provenance checks. Five official 2025-2026 research records were registered as motivation.
The full 170-script suite, memory 14, and heavy series stages were not launched. No critical endpoint
or thermodynamic observable moved.
The definitive archive contains 46 staged files plus manifest (47 validated members), is
1,575,644 bytes, and has SHA-256
`d167e544c7693b10d9526f66987da5a28cb71d610d5f3f4819831bf525f5a29e`; bundled documentation
omits archive self-references. The two earlier Wave-25 ZIPs are superseded and retained pending
user-approved deletion.
