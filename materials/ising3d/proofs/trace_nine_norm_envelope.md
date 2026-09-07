# Fixed-leading-form envelope for the open `2x4` trace-nine norm

Artifacts: `experiments/e252_trace_nine_norm_envelope.py`,
`results/spectral/trace_nine_lift_template.json.xz`,
`results/spectral/trace_nine_norm_envelope.json`, and the independent verifier
`tests/test_trace_nine_norm_envelope.py`. The e251 finite-projection theorem and
its frozen provenance remain in `proofs/trace_nine_projection.md`.

## 1. Result and boundary

**[THEOREM — exact norm envelope].** Work over `K=Q(q)` at values where the
normalized open-`2x4` trace targets are defined. The physical
`F4,...,F8` coefficient ideal has a monic 35-rule grevlex basis with the same
leading monomials as the homogeneous leading system. Its quotient has rank 96
and Hilbert vector

\[
(1,5,12,19,22,19,12,5,1).
\]

In that fixed basis, multiplication by `F9` is a `96x96` matrix of rational
functions of `q`. After clearing one common denominator per matrix row, let
`N(q)` be the resulting integer determinant numerator. The envelope proves

\[
\deg_q N\le 40693,
\qquad
\lVert N\rVert_1\le 2^{2343436}.                  \tag{1.1}
\]

The determinant denominator has degree at most `33017`. If `P` is the
primitive determinant numerator after polynomial cancellation, Gauss's lemma
makes `P` a primitive factor of `N`. The standard factor-height bound therefore
gives

\[
\max_j |[q^j]P|
\le 2^{\deg N}\lVert N\rVert_1
\le 2^{2384129}.                                  \tag{1.2}
\]

The q=2 modular Groebner basis and multiplication matrix are rebuilt from the
pinned trace rows and match the e251 certificate exactly: rank 96 and
determinant residue `1660951362` modulo `2147483647`. Hence the
characteristic-zero determinant is not identically zero.

**[RESOURCE OBSTRUCTION].** Dense scalar interpolation under the conservative
primitive-height envelope needs 40,694 nodes and at least 76,908 prime moduli
even under the optimistic ceiling of 31 information bits per safe prime. Thus
it needs at least 3,129,694,152 separate `96x96` determinant evaluations,
above the declared 10,000,000-evaluation gate. That interpolation was not run.

**[UNRESOLVED].** The primitive norm is not materialized. No real root is
isolated or classified, and open-`2x4` exceptional-set emptiness remains open.
The resource obstruction applies only to the named dense scalar-interpolation
plan; it does not rule out a fraction-free polynomial-matrix determinant,
structured reconstruction, or a sharper envelope. Nothing here is an
all-size, endpoint, or thermodynamic statement.

## 2. Homogeneous transformation certificate

Let `a=(a7,a6,a5,a4,a3)` be the five free coefficients after the fixed
`Q(0),Q(2),Q(3)` evaluations. Write

\[
F_i(a,q)=L_i(a)+R_i(a,q),\qquad i=4,\ldots,8,
\]

where `L_i` is the parameter-independent highest homogeneous part and
`deg_a R_i < deg_a L_i`. Their degrees are `(2,2,2,3,4)`. The XZ template
stores 35 exact rational identities

\[
G_j=\sum_{i=4}^8 H_{ji}L_i.                       \tag{2.1}
\]

Every `G_j` is monic and homogeneous. If `d_j=deg G_j`, every nonzero
multiplier `H_ji` is checked term by term to be homogeneous of degree
`d_j-deg L_i`; this graded check is what makes the lift below valid. The
stored basis is interreduced, its pure leading powers are `(2,2,3,5,9)`, and
literal standard-monomial enumeration gives 96 monomials with the stated
Hilbert vector. The producer validates every identity in (2.1) after reading
the compressed template; it does not trust stored leaders or dimensions.

The template is deliberately a deterministic certificate rather than a
second hand-maintained formula. `--build-template` derives the extended
Groebner transformations from the five leading forms, interreduces basis and
transformation rows together, writes canonical JSON through deterministic XZ,
and then reports both payload and compressed-file hashes.

## 3. Fixed-leading-form lift

Define

\[
\widetilde G_j
 =G_j+\sum_{i=4}^8 H_{ji}R_i
 =\sum_{i=4}^8 H_{ji}F_i.                         \tag{3.1}
\]

The graded condition above gives

\[
\deg_a(H_{ji}R_i)<d_j,
\]

so `\widetilde G_j` has the same monic leading term as `G_j`, with no division
by a q-dependent coefficient. Let `J=(\widetilde G_j)` and
`I=(F_4,\ldots,F_8)`. Equation (3.1) gives `J subseteq I`. The leading
monomials of `J` give a quotient of dimension at most 96. Conversely, the five
leading forms have no projective common zero, so the filtered/projective
Bezout argument from e248/e251 gives `dim_K K[a]/I=2*2*2*3*4=96`. The
surjection `K[a]/J -> K[a]/I` forces equality of the two dimensions and hence
`J=I`. Therefore the lifted rules are a Groebner basis of the physical ideal
with the fixed homogeneous initial ideal.

The implementation next orders the 35 lifts by ascending leader degree and
reduces each tail by earlier monic rules. These are subtraction and
multiplication operations only: no coefficient division occurs. It then
requires every final tail monomial to be one of the 96 standard monomials.
The resulting exact rational-function rules are bound by a canonical digest
in the JSON artifact.

## 4. Denominator, degree, and height propagation

For each exact rational coefficient, factor its primitive denominator into
integer q-polynomial atoms `A_k`. A bound record

\[
(d,h,e_0,e_1,\ldots)
\]

means that the coefficient has a representation

\[
\frac{P(q)}{\prod_k A_k(q)^{e_k}},\qquad
\deg P\le d,\qquad \lVert P\rVert_1\le 2^h.       \tag{4.1}
\]

The scalar denominator LCM is represented as a degree-zero atom, so (4.1)
always clears to an integer polynomial. Exact initialization computes the
actual numerator degree and coefficient `l1` norm. Multiplication adds degrees,
height exponents, and denominator valuations. Addition first takes the
componentwise maximum denominator valuation. If an operand is missing
`r_k` copies of `A_k`, its degree bound increases by
`r_k deg A_k` and its height exponent by
`r_k ceil(log2 ||A_k||_1)`; the two scaled numerators are then combined using
`max(h_left,h_right)+1`. These rules follow from
`||AB||_1 <= ||A||_1 ||B||_1` and
`||A+B||_1 <= ||A||_1+||B||_1`.

The producer propagates these records through the same ordered monic
reductions used for all 96 products `F9*b`. It records 8,557 bounded nonzero
matrix entries. For each matrix row it takes the componentwise maximum
valuation, scales every entry to that row denominator, and records the maximum
scaled numerator degree `d_r` and height exponent `h_r`. Leibniz then gives

\[
\deg N\le\sum_r d_r=40693,                        \tag{4.2}
\]

\[
\lVert N\rVert_1
 \le 96!\prod_r 2^{h_r}
 \le 2^{\lceil\log_2(96!)\rceil+\sum_r h_r}
 =2^{2343436}.                                    \tag{4.3}
\]


Polynomial cancellation is not coefficient-height monotone. Let `P` be the
primitive numerator after canceling factors shared by `N` and the product of
row denominators. Since `P` is a primitive factor of `N`, Mahler measure and
the standard coefficient bound for a degree-`m` factor give

\[
H(P)\le 2^m M(P)
     \le 2^{\deg N}M(N)
     \le 2^{40693}\lVert N\rVert_1
     \le 2^{2384129}.                              \tag{4.4}
\]

This additional factor is necessary: taking a polynomial factor can increase
coefficient height even when the dividend has smaller coefficients.

The sum of row-denominator degrees is `33017`. The artifact stores all 96 row
bounds, all 96 reduction-step counts, every denominator atom, and the final
rule digest; the displayed totals are recomputed rather than entered as free
metadata.

All q-polynomial denominator atoms come from the exact normalized targets.
The e251 proof shows their inherited physical denominators are strictly
positive for `q>1`; monic lifting and reduction introduce no new q-dependent
division. Row clearing may therefore add candidate factors algebraically but
does not discard a physical point.

## 5. Certified reconstruction stop

A dense polynomial of degree at most 40,693 needs 40,694 values per prime.
For signed CRT uniqueness under the conservative primitive bound (4.4), the
combined modulus must be strictly larger than `2^(2384129+1)`. The e251
arithmetic admits primes below `2^31`; even the optimistic information ceiling
of 31 bits per prime implies

\[
\left\lfloor\frac{2384129+1}{31}\right\rfloor+1
 =76908
\]

prime moduli are necessary. Multiplying by the node count yields the certified
lower bound 3,129,694,152 determinant evaluations. Because this already
exceeds the explicit cap, the producer records
`NOT_RUN_CERTIFIED_CRT_LOWER_BOUND_EXCEEDS_CAP` and stops. Apparent degree
stabilization at a few nodes or primes would not certify reconstruction and is
not used.

The highest-value continuation on this same gate is a fraction-free or modular
polynomial-matrix determinant exploiting the fixed 96-dimensional template,
or a rigorously sharper denominator/height structure. Only after an exact
primitive norm exists can every `q>1` root and its repeated-root-safe shifted
`H0/H1` coefficient branch be classified.

## 6. Mutation and independent verification

Six in-memory mutations exercise distinct load-bearing claims:

1. change one transformation coefficient, breaking (2.1);
2. replace the pure leader `a3^9` by `a3^10`, changing the standard-monomial
   count;
3. lower the stored determinant-degree envelope by one and compare it with the
   independently recomputed row sum;
4. lower the stored primitive-height envelope by one and compare it with the
   independently recomputed factor-height bound;
5. replace the inherited nonzero q=2 norm witness by zero;
6. replace the inherited `q^(-20)` trace-four normalization by `q^(+20)`, whose
   q=2 ratio is `2^40`.

All six are rejected and their witnesses are stored. The verifier imports no
e252 producer code. It rebuilds the physical incidence through the independent
e251 verifier, checks every XZ transformation and graded identity, uses a
separate highest-grevlex-term scan rather than the producer's heap reducer,
requires every defining `F4,...,F8` equation to reduce to zero under the final
lift, and reconstructs all exact lifts and bound records. Independently of
that replay, it recomputes the determinant totals from the artifact's 96 row
bounds, validates each stored denominator atom directly, and rebuilds both
e251 physical witnesses at `q=2` and `q=5/3`. It then compares the complete
check and data payloads, validates every current source hash, and replays the
six mutation witnesses.

From the `math/` repository root, the deterministic commands are:

```sh
nice -n 19 ./.venv/bin/python -I -B ising3d/experiments/e252_trace_nine_norm_envelope.py --build-template
nice -n 19 ./.venv/bin/python -I -B ising3d/experiments/e252_trace_nine_norm_envelope.py
nice -n 19 ./.venv/bin/python -I -B ising3d/tests/test_trace_nine_norm_envelope.py
```

The actual interpreter is `math/.venv/bin/python` at Python `3.14.3`. Its
dependency contract is `math/requirements-freeze.txt`, including
`python-flint==0.9.0` backed by FLINT `3.6.0`; the producer and verifier now
compare every frozen distribution against the installed environment and pin
the freeze alongside the inherited `ising3d/uv.lock`.
