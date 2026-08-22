# Layer-group spectral route: an exact varying-conjugator obstruction

Artifacts: `experiments/e188_gaussian_group.py`,
`experiments/e189_varying_conjugator.py`,
`experiments/e190_layer_group_spectral.py`, and
`results/spectral/layer_group_spectral.json`.  Independent verifier:
`tests/test_layer_group_spectral.py`.

Reproduce from the repository root:

```bash
.venv/bin/python experiments/e190_layer_group_spectral.py
.venv/bin/python tests/test_layer_group_spectral.py
```

Every mathematical statement below is tagged **[THEOREM]**, **[LEMMA]**,
**[COMPUTATION]**, **[EXTERNAL]**, **[CONJECTURE]**, or **[UNRESOLVED]**.

## 0. Outcome and exact scope

**[THEOREM — exact obstruction].**  The proposed rigidity bridge is false:
pointwise membership

```text
exp(aA) exp(bB) in union_g g (C^times rho(Spin(2n,C))) g^(-1)
```

on an open set of parameters does **not** force `Lie(A,B)` into one fixed
conjugate of the quadratic algebra.  Section 4 gives a one-mode spin-group
witness for which every point in a nonempty complex-open two-parameter set is
conjugate into `Spin(2,C)`, while the two infinitesimal generators span
`sl(2,C)`.  Section 5 gives the exact varying-conjugator velocity that absorbs
the entire transverse derivative at `b=0`.

**[UNRESOLVED — requested spectral theorem].**  This front does **not** prove
that every branching layer has a non-Gaussian spectrum at all but finitely many
physical isotropic couplings.  No all-size exceptional set is obtained:
`E_G` remains **undefined/unresolved**, not empty.  This is a proof-method
obstruction, not a Gaussian branching-layer example and not a solution of the
three-dimensional Ising model.

**[COMPUTATION — controls only].**  At `t=1/3`, one-prime exact modular
certificates independently reproduce the established `2x3` and `2x4` spectral
no-gos: respectively `r_Q>=1631>665` and `r_Q>=23311>6305`.  These are finite
controls for the machinery, not the primary method and not an all-size
induction.

## 1. The two target formulations

**[LEMMA — definitions].**  Let

```text
F_n = Lambda^bullet C^n,                         dim F_n = 2^n,
q_n = span_C { gamma_i gamma_j : 1<=i<j<=2n }  ~= so(2n,C),
H_n = rho(Spin(2n,C)) <= GL(F_n),
Htilde_n = C^times H_n.
```

Here `rho` is the full spinor/Fock representation; it contains both half-spin
parity sectors.  The scalar extension `Htilde_n` is required because the
physical transfer operator is defined only up to a nonzero global scalar.  Say
that an invertible operator is *group-Gaussian in some basis* when it belongs
to

```text
Sat(Htilde_n) = union_{C in GL(F_n)} C Htilde_n C^(-1).             (1)
```

**[LEMMA — spectral target].**  A *nonzero full `n`-mode subset-product
multiset* is the `2^n`-slot multiset

```text
Sigma_n(a;u_1,...,u_n)
  = { a product_{i in S} u_i : S subseteq {1,...,n} },
a,u_i in C^times,                                                   (2)
```

with one slot for every subset, including multiplicity.  Requiring the
parameters to be nonzero matches the invertible physical transfer operator and
the group in (1).  Allowing a zero `a` or `u_i` defines a larger spectral
closure, not a subset of a matrix group.

## 2. Exact implication directions

**[EXTERNAL — standard spin-torus facts].**  A semisimple element of a complex
spin group is conjugate into a maximal torus.  On the full spinor module the
torus weights can be written

```text
c product_{i=1}^n r_i^(epsilon_i),
epsilon_i in {-1,+1},  r_i in C^times.                              (3)
```

For a general group element, multiplicative Jordan decomposition and
simultaneous triangularization show that its eigenvalues are those of its
semisimple part; the commuting unipotent part changes Jordan blocks but not
eigenvalues.

**[LEMMA — forward direction, no converse silently used].**  If `W` is
conjugate into `Htilde_n`, then `spec(W)` is a nonzero full `n`-mode
subset-product multiset.  No diagonalizability assumption on `W` is required
for this direction.

**[LEMMA — proof].**  Apply (3) to the semisimple part and put

```text
a = c product_i r_i^(-1),        u_i = r_i^2.
```

The map `epsilon_i -> (epsilon_i+1)/2` bijects the `2^n` spin weights with the
`2^n` subsets, and (3) becomes exactly (2), slot by slot. `[]`

**[THEOREM — converse with its necessary hypotheses].**  If `W` is
**invertible and diagonalizable** and `spec(W)=Sigma_n(a;u_1,...,u_n)` with all
parameters nonzero, then `W` is conjugate into `Htilde_n`.

**[THEOREM — proof].**  Choose square roots `r_i^2=u_i` and put
`c=a product_i r_i`.  The torus element with weights (3) has exactly the
multiset (2).  Two diagonalizable complex matrices with the same eigenvalue
multiset, including algebraic multiplicities, are similar.  Therefore `W` is
similar to this scalar spin-torus element. `[]`

**[THEOREM — why the unrestricted spectral converse is false].**  Spectrum
alone does not determine Jordan structure.  The nontrivial Jordan block

```text
J_2(1) = [[1,1],[0,1]]
```

has the one-mode subset-product multiset `{1,1}` (`a=u_1=1`), but it is not
conjugate into the scalar `Spin(2,C)` torus, whose matrices are diagonalizable.
Likewise, a subset product with a zero parameter cannot be conjugate into any
matrix group.  Thus the forward implication is unconditional, while the
reverse implication needs precisely the invertibility and diagonalizability
used above.

**[LEMMA — physical operators meet the converse hypotheses].**  For real
`a,b`, with Hermitian Pauli sums `A` and `B`,

```text
V(a,b) = exp(aA) exp(bB)
```

is similar to

```text
S(a,b) = exp(aA/2) exp(bB) exp(aA/2)
```

via `V=exp(aA/2) S exp(-aA/2)`.  The sandwich `S` is Hermitian positive
definite, hence invertible and diagonalizable with positive real spectrum.

**[THEOREM — exact relation for the physical target].**  Consequently, for the
real physical layer and up to its nonzero scalar normalization,

```text
V(a,b) is group-Gaussian in some basis
    <=> spec(V(a,b)) is a full nonzero n-mode subset-product multiset.          (4)
```

The left-to-right implication holds more generally; the right-to-left
implication in (4) is valid here because of the preceding positivity lemma.

## 3. The physical family is a proper curve

**[LEMMA — exact curve equation].**  Put

```text
t = tanh(a/2),          x = exp(a),          y = exp(2b).
```

The physical isotropic duality relation is

```text
x = (1+t)/(1-t),        y = (1+t^2)/(2t).
```

Eliminating `t=(x-1)/(x+1)` gives the proper curve

```text
F_phys(x,y) = (y-1)x^2 - (y+1) = 0.                              (5)
```

The producer verifies (5) in `Q(t)` and verifies that `F_phys` has the four
nonzero coefficients `x^2 y - x^2 - y - 1`; hence it is not the zero ambient
polynomial.

**[THEOREM — ambient noncontainment is insufficient].**  Even if one proved
that the full two-parameter image `(x,y) -> V(x,y)` was not contained in
`Sat(Htilde_n)`, that statement alone would not show that the restriction (5)
was not contained.  A nonzero ambient obstruction polynomial can vanish
identically on the physical curve: `F_phys` itself is an exact witness, and so
is every nonzero multiple of it.

**[LEMMA — what a finite exceptional set would require].**  A polynomial
spectral obstruction `P_G(x,y)` yields an all-but-finitely-many theorem along
the physical curve only after proving

```text
P_G(x(t),y(t)) != 0 in Q(t),
```

or equivalently that the relevant class of `P_G` is nonzero modulo
`F_phys`.  Its nonzero cleared univariate numerator then has finitely many
roots.  Ambient two-parameter noncontainment supplies no such nondivisibility.
This is separate from, and in addition to, the varying-conjugator failure below.

## 4. Counterexample to the proposed rigidity theorem

**[THEOREM — exact one-mode witness].**  On `C^2` define

```text
A_0 = H = [[1,0],[0,-1]],          B_0 = X = [[0,1],[1,0]].
```

The one-mode spin image is

```text
T = rho(Spin(2,C)) = { diag(r,r^(-1)) : r in C^times }.
```

Its scalar extension is the group of invertible diagonal matrices, and its
`GL(2,C)` conjugacy saturation is the set of invertible diagonalizable `2x2`
matrices.

**[LEMMA — open family in varying conjugates].**  In multiplicative parameters
`r=exp(a)`, `u=exp(b)`, write

```text
c=(u+u^(-1))/2,       s=(u-u^(-1))/2,
M(r,u)=exp(aH)exp(bX)
      = [[r c,       r s],
         [r^(-1)s, r^(-1)c]].                                      (6)
```

Then

```text
det M = 1,
tr M = (r+r^(-1))(u+u^(-1))/2,
Delta(r,u) = (tr M)^2-4.                                           (7)
```

On the nonempty complex-open set

```text
Omega = { (r,u) in (C^times)^2 : Delta(r,u) != 0 },                 (8)
```

`M` has two distinct nonzero eigenvalues `lambda,lambda^(-1)`, so it is
conjugate to `diag(lambda,lambda^(-1)) in T`.  The exact non-diagonal control
`(r,u)=(2,2)` gives

```text
M = [[5/2,3/2],[3/8,5/8]],       Delta=369/64 != 0.
```

The open set also contains regular points with `b=0`, for example `(r,u)=(2,1)`.

**[LEMMA — real-open strengthening].**  For real additive parameters,

```text
tr M = 2 cosh(a) cosh(b) >= 2.
```

Equality occurs only at `(a,b)=(0,0)`, where `M=I`; everywhere else the two
eigenvalues are distinct.  Thus every real point of the two-parameter family is
pointwise conjugate into `T`, including the identity.

**[LEMMA — the pair is not jointly quadratic].**  Exact brackets give

```text
[H,X] = [[0,2],[-2,0]],
[H,[H,X]] = 4X,
[X,[H,X]] = -4H.
```

Therefore `Lie(H,X)=sl(2,C)` has dimension three.  It cannot be conjugate into
the one-dimensional `Lie(T)`; even adding scalar matrices gives only dimension
two.

**[THEOREM — refutation].**  Equations (6)--(8) give a two-parameter family
lying pointwise in varying conjugates of the spin subgroup on a nonempty open
set, while its generator pair is not jointly conjugate into quadratics.  Hence
the proposed implication in §0 is false for the spin subgroup itself, not merely
for an unrelated abstract subgroup. `[]`

## 5. Why differentiation at `b=0` does not repair it

**[LEMMA — tangent space of a conjugacy saturation].**  For an algebraic
subgroup `H<=G`, differentiate

```text
Phi(g,h)=g h g^(-1).
```

At `(1,h)`,

```text
d Phi_(1,h)(Y,Z) = [Y,h] + Z,        Z in T_h H.                    (9)
```

Thus motion of the conjugator contributes the entire conjugacy-orbit tangent
`[Lie(G),h]`; differentiation does not merely return `T_h H`.

**[LEMMA — exact absorption in the witness].**  At `b=0`, put

```text
D=diag(r,r^(-1)),
partial_b M|_(b=0) = D X,
Y(r) = 1/(r^2-1) [[0,-r^2],[1,0]].                                 (10)
```

A direct multiplication gives

```text
[Y(r),D] = D X,             r^2 != 1.                              (11)
```

Thus the entire `B_0` derivative is a conjugacy-orbit velocity in (9), with zero
torus velocity.  At the regular point `r=2`,

```text
Y(2)=[[0,-4/3],[1/3,0]],
```

so the conjugator varies at a finite rate and exactly absorbs the apparently
nonquadratic direction.

**[THEOREM — the identity is precisely singular].**  The velocity `Y(r)` has a
pole at `r^2=1`.  Therefore a proof that assumes one bounded differentiable
choice of diagonalizing conjugator through the degenerate identity has assumed
the missing conclusion.  Pointwise conjugators exist at the identity, but no
such regular choice is forced.  The phrase “the conjugator cannot vary fast
enough” fails here in the strongest possible way: the required rate is finite
at every regular `b=0` point and becomes unbounded at the degeneracy where the
naive tangent argument would try to anchor it.

## 6. Consequence for the layer Lie algebra

**[THEOREM — fixed-basis scope].**  The exact local-term dimensions in
`proofs/clifford_grade_classification.md` and any lower bound on
`Lie(A,B)` can rule out a **single fixed** matrix `C` satisfying

```text
C^(-1) A C, C^(-1) B C in q_n.
```

They cannot, by dimension alone, rule out separate conjugators `C(a,b)` for the
individual matrices `exp(aA)exp(bB)`.  The witness has exactly this separation:
its generated algebra is too large for one fixed torus, although every regular
family member is individually torus-conjugate.

**[LEMMA — `ad_A` decomposition does not by itself separate grades].**  For
the layer, `A=sum_v X_v` has Clifford grade two, so `ad_A` preserves each
Clifford-grade subspace.  Iterated brackets and a Vandermonde argument can
separate distinct `ad_A` eigenvalues occurring in `B`, but preservation alone
does **not** imply that different Clifford grades have disjoint eigenvalue
sets.  Thus this observation neither proves that `B_2` and `B_high` can be
separated nor supplies any equality or regularity relation among the pointwise
spectral conjugators `C(a,b)`.  Even a successful fixed-algebra lower bound
would not contradict (9)--(11) and would not prove a spectral no-go.

**[THEOREM — independence from the contemporaneous dimension front].**  The
obstruction proved here is logical rather than numerical: regardless of how
large a correct two-generator dimension lower bound is, it excludes only joint
fixed-basis quadraticity until an additional theorem controls the varying
conjugator.  Accordingly, no unverified size threshold from the concurrent
two-generator front is imported into this proof.

## 7. Exact finite controls

**[COMPUTATION — raw construction].**  For each open grid below, the producer
builds the exact integer multiple of

```text
R = P_t D_q P_t,       t=1/3, q=5/3,
```

then reconstructs the characteristic polynomial and the unordered-slot pair
polynomial over `F_1000003` from exact power traces, Newton identities, and the
Cayley--Hamilton trace recurrence.  It computes the monic derivative gcd by
Euclid and checks that the gcd divides both inputs.  Peak RSS is about `103 MB`
and producer process time is below `7 s` on the recorded host.

**[LEMMA — modular direction].**  The characteristic-zero pair polynomial is
monic.  At a prime of good reduction,

```text
deg gcd_Q(C_2,C_2') <= deg gcd_Fp(C_2 mod p,C_2' mod p),
```

so a modular degree `g_p` proves only

```text
r_Q >= binom(2^n,2)-g_p.                                           (12)
```

This one-sided lower-bound direction is the only direction used.

**[COMPUTATION — control table].**

| layer | `n` | edges | pair slots | `g_p`, `p=1000003` | proved `r_Q>=` | Gaussian ceiling `3^n-2^n` |
|---|---:|---:|---:|---:|---:|---:|
| open `2x3` | 6 | 7 | 2016 | 385 | **1631** | 665 |
| open `2x4` | 8 | 10 | 32640 | 9329 | **23311** | 6305 |

**[COMPUTATION — independent verification].**  The standalone test imports no
producer.  It constructs `P_t` entry by entry from Hamming distances, constructs
`D_q` from raw grid edges, multiplies `P_t D_q P_t` directly in `F_1000003`, and
independently repeats the power-trace/Newton/Euclid calculation.  It reproduces
both gcd degrees and both one-sided bounds, verifies symmetry and raw edge
counts, stays below `51 MB` peak RSS, and prints `PASS test_layer_group_spectral`.
The controls agree with the prior certificates in
`proofs/allsize_gaussian.md`; no prior producer conclusion is imported by this
test.

## 8. Final statement and next viable routes

**[THEOREM — route closed].**  The combination

```text
large Lie(A,B)
+ differentiation of pointwise Gaussian conjugacies
=> fixed quadratic conjugate
```

is invalid.  The exact failure term is `[Y,V]` from the derivative of the
conjugator, and equations (10)--(11) show it can absorb the whole transverse
direction.  Separately, ambient two-parameter noncontainment does not control
the proper physical curve (5).

**[UNRESOLVED — exceptional set].**  No polynomial `rho_G(t)` and no finite
exceptional set are produced for arbitrary branching layers.  The two named
control layers remain certified non-Gaussian at `t=1/3` (and have their prior
named-graph generic consequences), but this front establishes no all-size
quantifier.

**[UNRESOLVED — what a successful continuation must add].**  A group/spectral
continuation must provide at least one ingredient absent here:

1. a conjugacy-invariant spectral polynomial `P_G(x,y)` proved nonzero modulo
   `F_phys` for every branching layer; or
2. additional hypotheses, actually satisfied by the Ising family, that force
   the pointwise conjugators to come from one parameter-independent conjugator.

A Lie-dimension bound without one of these ingredients cannot decide the
physical subset-product spectrum.
