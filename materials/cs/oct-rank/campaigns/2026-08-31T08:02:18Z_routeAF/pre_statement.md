# Route AF pre-statement

Timestamp: 2026-08-31T08:02:18Z. Owner: agent `OctRankRouteAF`.
This file is the first file in this campaign directory and is committed before
campaign computation. A prior session contained unrecorded exploratory
numerical design notes; no result from that exploration is adopted here. Every
reported result must be regenerated after this pre-statement commit.

## Fixed convention and targets

Use Cayley--Dickson octonions in the ordered basis
`(1,i,j,k,l,il,jl,kl) = (e0,...,e7)`, with
`(a,b)(c,d)=(ac-conj(d)b, da+b conj(c))`.  For `x`, the left-multiplication
matrix has column `b` equal to the coordinates of `x*e_b`, so
`L_x[c,b]=(x*e_b)_c`.  The Route-F tensor is the `3 x 8 x 8` tensor whose
three matrix slices are `(L_1,L_i,L_j)` in this convention.

**Job F.** Decide whether this tensor has real tensor rank 13 or 14.  The
inherited lower chain is 13 and the blockwise `tau+tau` upper bound is 14.
Success downward requires an exact-rational decomposition or a certified
interval-Newton/Krawczyk existence certificate for 13 real rank-one terms,
plus replay of the inherited lower chain.  Success upward requires a valid
rank argument certifying rank at least 14; a failed numerical search is not a
lower bound.

**Job A.** Audit, from primary sources, the precise hypotheses and conclusion
of the proposed Strassen/Blaser three-slice commutator inequality.  Separately
seek an algebraic rank proof for
`[L_conj(u)L_v,L_conj(u)L_w]` on every linearly independent real octonion
triple.  A finite sweep is evidence only.  The mandatory `n=4` control is the
quaternion tensor `tau=(L_1,L_i,L_j)`, whose exact rank is 7; any proposed
formula returning 8 for this control is refuted or missing a hypothesis and
must not be used.

## Anchors and adjudication

Before trusting either instrument, replay these known values at startup:

1. complex multiplication rank 3: exact three-term Gauss witness and exact
   lower-bound anchor;
2. quaternion multiplication rank 8: exact eight-term witness and inherited
   certified lower chain;
3. `tau=(L_1,L_i,L_j)` rank 7: the frozen Krawczyk certificate and lower chain;
4. the Route-F tensor's exact block identities and the certified 14-term
   upper construction.

For Job F, plant rank-13-feasible synthetic controls with the same tensor
shape before interpreting any zero/failure count.  Numerical residuals,
least-squares convergence, and float singular values are
`COMPUTATIONAL-EVIDENCE` only.  `MACHINE-VERIFIED` requires exact rational
identities or certified interval bounds with positive width/margins end to
end.  The stored decimal factors of a Krawczyk certificate are approximate;
they are never described as an entrywise-exact witness.

## Pre-registered computation budget

Job F numerical discovery may use at most 32 deterministic seeds for a
rank-13 nonlinear least-squares search, at most 5,000 accepted optimization
steps per seed, and a success trigger of relative Frobenius residual at most
`1e-11` with bounded factors.  Any candidate crossing that trigger is
escalated to exact/certified replay; the campaign stops discovery once a valid
certificate is obtained.  Exact certification may use rational reconstruction,
Jacobian rank checks, projection to a square nonsingular subsystem, and
Krawczyk/interval Newton at any rational precision needed.  Failure to certify
is reported with the best residual, factor norms, seeds, and the exact failed
margin; it makes no rank lower-bound claim.

Job A may use exact symbolic identities and determinant/rank calculations over
`fmpq`; any universal claim must end in a proof identity or a certified
nonvanishing formula on the exact declared domain.  Primary-source audit must
record edition/version, page or line ranges, and verbatim relevant hypotheses.

## Pre-registered outcomes

- **F13-certified:** exact/certified rank-13 upper plus lower 13 replay, hence
  Route-F tensor rank exactly 13.
- **F14-certified:** algebraic rank-at-least-14 proof plus existing upper 14,
  hence exact rank 14.
- **FAILURE TO CERTIFY:** neither side closes; report quantified numerical and
  certificate obstructions without changing the `[13,14]` interval.
- **A-live:** provenance verified, tau control passed, and a universal rank
  argument discharges the commutator hypothesis.
- **A-dead/unverified:** source theorem does not imply the proposed formula,
  tau control refutes the reading, or the universal rank argument is absent.

A purported change to the published `18 <= R_R(T_O) <= 25` window is escalated
to Main before recording.

## Rule-7 prospective scope

The campaign sweeps exactly the fixed Route-F tensor `(L_1,L_i,L_j)` in the
8-dimensional Cayley--Dickson convention above, rank target 13, 32 seeds,
5,000 accepted steps per seed, and whatever certified interval radius is
reported for the selected candidate.  It does not search other octonion
triples, other ranks of the full `8 x 8 x 8` multiplication tensor, complex
rather than real decompositions, border rank, non-CP constructions, or lower
bounds by finite-field/rational enumeration.  Job A concerns all linearly
independent triples `(u,v,w) in R^8` only if a universal proof is found;
otherwise its exact scope is the explicitly listed controls and any finite
sweep, with no extrapolation.
