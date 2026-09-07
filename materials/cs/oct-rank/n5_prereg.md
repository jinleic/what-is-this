# N5 pre-statement — gate `n5-direct-strassen-q`: exact complex rank of quaternion multiplication

Campaign: mint only after this file is committed, with
`python3 ../../scripts/campaign.py init --gate n5-direct-strassen-q --prereg n5_prereg.md --agent OctRankN5`
from `cs/oct-rank/`. Date: 2026-09-04. No claim-relevant calculation may
run before that init.

## 0. Invalid-prereg history and clean cutover

Runs `20260904T042512Z_0b95f53e_d15e6b73a0ae` and
`20260904T043034Z_334d11f6_e28d8f892300` are preserved
FROZEN-INCONCLUSIVE invalid-prereg history. Neither produced an N5 witness or
adjudication. This remint removes both sources of ambiguity in the second
attempt: there is one output orientation, and no slice-algebra/S0 pullback is
part of witness construction.

## 1. Registered claim and evidence labels

The registered claim is

`rank_C(T_H) = 7`,

where `T_H` is the 4 x 4 x 4 structure tensor of quaternion multiplication in
the ordered basis `(1,i,j,k=ij)`. Every equality restatement must carry the
composite label

**upper bound MACHINE-VERIFIED / lower bound CITED-DEPENDENCY
(source-locked)**.

* Upper bound: an explicit seven-term CP decomposition over Q(i), constructed
  in the run and checked exactly against all 64 entries of `T_H`. This part is
  MACHINE-VERIFIED.
* Lower bound: the Alder-Strassen inequality `L(A) >= 2 dim(A) - t(A)`, with
  `t(A)` the number of maximal two-sided ideals, applied to the simple
  four-dimensional C-algebra `M_2(C)`, hence `7`. This theorem and its
  connection to multiplication-tensor rank are consumed as a
  CITED-DEPENDENCY, not rederived or relabeled machine-certified.

Source lock: A. Alder and V. Strassen, “On the algorithmic complexity of
associative algebras,” *Theoretical Computer Science* 15 (1981), DOI
`10.1016/0304-3975(81)90070-0`. The pre-existing source-lock records are bound
by SHA-256:

* `scratch/n5_as_source_lock/crossref_record.md`:
  `d66a28dadea33518a1effa08d596bd53e1e1c9b6630e884d0063fed577787be8`;
* `scratch/n5_as_source_lock/source_lock.md`:
  `9e1e88396494c422b186df05ebc163ba4e28ec390c2a5d3ef4e1df5d57bf7849`.

The full paper was not accessible in the prior source-lock session; that
limitation is why the lower bound remains CITED-DEPENDENCY.

This gate makes no new real-rank claim. In particular, the real
`{13,14}` frontier does not move, and `18 <= R_R(T_O) <= 25` is untouched.
“Absence of a 13-witness is NOT evidence for 14; absence of an impossibility
argument is NOT evidence for 13.”

## 2. Fixed direct algebra construction

All indices below are zero-based in code. Human matrix-unit subscripts in the
simplicity proof are one-based.

### 2.1 Quaternion tensor

Construct the Hamilton table directly in `(1,i,j,k)` and set, without any
slice transpose or change of action,

`mult_Q[p,q,r] = T_H[p][q][r]`.

The instrument asserts the 16-nonzero table anchor. `mult_Q` is the sole
target tensor used in witness construction and final substitution.

### 2.2 Exact isomorphism to M2(Q(i))

Let the scalar field imaginary unit be denoted `ii`, and fix

```
I2 = [[1,0],[0,1]]
Xi = [[0,1],[-1,0]]
Xj = [[0,-ii],[-ii,0]]
```

with ordered quaternion image basis

`B = [I2, Xi, Xj, XiXj]`.

The instrument verifies exactly that `Xi^2=Xj^2=-I2` and
`Xi Xj = -Xj Xi`. Flatten matrices in row-major E-order
`(E11,E12,E21,E22)`. Let `C[a,p]` be the E-coordinate of `B[p]`; thus the
columns of `C` are the four listed basis matrices, and let
`Cinv = C^-1`.

Let `mult_std[a,b,d]` be standard 2 x 2 matrix multiplication in that
row-major E-basis. The exact isomorphism identity is fixed as

```
mult_Q[p,q,r]
  = sum_{a,b,d} C[a,p] C[b,q] mult_std[a,b,d] Cinv[r,d].
```

The orientation is **`Cinv[r,d]`**, because an E-coordinate output is
converted back to B-coordinates by left multiplication with `C^-1`. All 64
identities, `C Cinv = Cinv C = I4`, and the defining quaternion relations are
asserted before witness assembly. There is no fallback orientation.

### 2.3 Strassen terms and exact output solve

Fix only the seven classical input linear-form pairs, in row-major E-order:

```
((1,0,0,1),  (1,0,0,1))
((0,0,1,1),  (1,0,0,0))
((1,0,0,0),  (0,1,0,-1))
((0,0,0,1),  (-1,0,1,0))
((1,1,0,0),  (0,0,0,1))
((-1,0,1,0), (1,1,0,0))
((0,1,0,-1), (0,0,1,1))
```

No output vector is transcribed. For each output coordinate `d`, the run forms
the 16 x 7 exact coefficient matrix
`D[(a,b),k] = u_k[a] v_k[b]` and solves

`D w[:,d] = mult_std[:,:,d]`

by deterministic exact Gaussian elimination over Q (embedded in Q(i)). It
asserts consistency, rank seven, uniqueness, and exact 64/64 reconstruction
of `mult_std`.

Transport each solved term by exactly

`(u'_k, v'_k, w'_k) = (C^T u_k, C^T v_k, C^-1 w_k)`.

Equivalently, `w'_k[r] = sum_d Cinv[r,d] w_k[d]`. The run then substitutes
these seven transported terms **directly** into `mult_Q=T_H` and requires
64/64 exact matches and zero mismatches. No `S0`, cyclic vector, slice
algebra, conjugated frame, opposite algebra, or asymmetric matrix-factor
pullback enters this construction.

## 3. Simplicity and cited lower bound

The run records the following correct algebra proof supporting the cited
theorem’s `t=1` input. If a nonzero two-sided ideal `I` of `M_2(C)` contains
`M != 0`, choose `M_ab != 0`. Then

`M_ab^-1 E_1a M E_b2 = E_12`.

Hence `E_12` lies in `I`, and for every `i,j`,

`E_ij = E_i1 E_12 E_2j`

lies in `I`. Thus `I=M_2(C)`; `(0)` is the only proper and hence the unique
maximal two-sided ideal, so `t=1`. The instrument spot-verifies the relevant
matrix-unit identities for every `a,b,i,j` over the exact backend. The
Alder-Strassen dependency then supplies `2*4-1=7`; the run labels this
CITED-DEPENDENCY throughout.

Any optional slice-algebra calculation is forbidden from witness construction
and adjudication. The direct `mult_Q` isomorphism above is sufficient.

## 4. Registered controls and independent audit

1. **A7 accept:** the newly transported seven terms must reproduce exactly all
   64 entries of `T_H`; term count is exactly seven and no factor vector is
   zero.
2. **R7 reject:** increment only `T_H[0][0][0]` by one. The same seven terms
   must fail, with the exact mismatch count and coordinates recorded.
3. **A8 anchor:** read the frozen N4 artifact
   `campaigns/20260904T035341Z_cc4a6ac5_e25c6371c8ce/witness_tower_qi.json`,
   require SHA-256
   `70dfc93c6496a0390d1a9a86af668e10a28e6b65dbe8e5d617fba2964272fc5b`,
   and independently accept its `TH_8term_orig` value on all 64 entries.
   The corrupted target must reject that anchor too. This is a control only;
   it does not enter the seven-term witness.
4. **Independent re-substitution:** after the authoritative instrument emits a
   serialized seven-term witness, a separate stdlib-only verifier that does
   not import the instrument rebuilds the Hamilton table from an independent
   basis-index/sign table, parses every rational complex coefficient, checks
   all 64 entries, repeats the single-coordinate corruption rejection, and
   writes a separate audit JSON containing the witness SHA-256.

A control false accept, false reject, dependency-hash mismatch, dimension
mismatch, nonunique solve, or convention assertion failure yields
FROZEN-INCONCLUSIVE (INVALID INSTRUMENT/CONSTRUCTION), never certification.

## 5. Runtime and artifact discipline

* Python stdlib only; exact `fractions.Fraction` pairs represent Q(i); no
  floats in any claim-relevant arithmetic or serialized coefficient.
* Launch with `PYTHONDONTWRITEBYTECODE=1 nice -n 10`; set
  `sys.dont_write_bytecode=True` before non-stdlib local imports (there are
  none); assert the environment flag and record process priority.
* Set `RLIMIT_CPU` soft and hard to an effective cap no larger than 7200 CPU
  seconds. Check the budget at every stage boundary. On an orderly budget
  breach, persist a budget-abort JSON for frozen inconclusive closeout.
* One deterministic pass, no parameter adjustment, fallback, or search.
* Artifacts: byte-identical `pre_statement.md`, provenance with prereg commit
  and source/dependency hashes, `n5_instrument.py`, `n5_results.json`,
  `n5_witness_qi.json`, independent verifier and audit, `VERDICT.md`, and the
  freeze-generated `sha256s.txt` ledger. Preserve every failure artifact.

## 6. Verdict rule and scope

**FROZEN-CERTIFIED** iff every construction assertion and registered control
passes, authoritative direct substitution reports 64/64, the independent
verifier reports 64/64 and rejects the corruption, all hashes match, and the
source dependency is labeled correctly. The certified statement is only:

`rank_C(T_H)=7` — **UB MACHINE-VERIFIED / LB CITED-DEPENDENCY
(source-locked)**.

Otherwise close FROZEN-INCONCLUSIVE with the exact defect or resource reason.

Covered: this one complex-rank equality, its direct seven-term upper witness,
the exact direct isomorphism checks, matrix-unit simplicity support, the
source-locked cited lower bound, and the named controls. Not covered: any real
rank, the real `{13,14}` frontier, `T_O`, border rank, rederivation of the
Alder-Strassen theorem, or a new statement about `R_R(T_O)`.
