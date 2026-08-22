# Pair-product collision obstruction for the exact `2x3` layer spectrum

Artifacts: `experiments/e107_pair_product_obstruction.py`,
`tests/test_pair_product_obstruction.py`, and
`results/spectral/pair_product_obstruction.json`.

**Status tags.** **[THEOREM]** and **[LEMMA]** are proved below. **[COMPUTATION]**
means an exact finite-field/integer computation recorded in the artifact and
independently rebuilt by the standalone test. **[UNRESOLVED]** marks limits of
this certificate.

This is an independent, **ordering-free** certificate for the exact-coupling
statement of Theorem S. It never sorts eigenvalues, invokes inertia, selects the
lowest eigenvalues, or assumes simplicity.

## 1. Statement and scope

**[THEOREM P — pair-product no-go].** At

\[
t=\tanh(K^*/2)=\frac13,\qquad e^{2K}=\frac53,
\]

the transfer operator of the open `2x3` Ising layer is **not** a full
six-mode Gaussian subset-product spectrum

\[
\left\{a\prod_{i=1}^6u_i^{\epsilon_i}:
  \epsilon\in\{0,1\}^6\right\}
\tag{1}
\]

as a multiset of 64 eigenvalues. The conclusion permits arbitrary complex
parameters in (1), hence in particular rules out the positive-real parameters
of a positive-spectrum fermionic Gaussian.

The e38 construction gives a rational symmetric positive-definite matrix `R`
which differs from the physical transfer operator only by similarity and a
nonzero global scalar. Both operations preserve pair-product collision
multiplicities; a scalar merely replaces `a` in (1). It therefore suffices to
work with `R`.

The theorem is scoped only to a **full** six-mode spectrum at this one rational
coupling and this one layer. The exclusions are collected in sec. 7.

## 2. Pair products use distinct *occurrences*, not distinct values

Let \(\lambda_1,\ldots,\lambda_{64}\) be any listing of the eigenvalues of
`R` in an algebraic closure, repeated by algebraic multiplicity. Define

\[
 C_2(z)=\prod_{1\le i<j\le64}(z-\lambda_i\lambda_j).
\tag{2}
\]

Thus `i < j` means two distinct **eigenvalue slots**. If, for example,
\(\lambda_i=\lambda_j\) numerically, the pair remains included. There are
\(\binom{64}{2}=2016\) factors. This convention is essential: it handles all
spectral degeneracies without choosing an ordering of the distinct values.

### [LEMMA 1 — Gaussian pair-product collision bound]

If the eigenvalue multiset has the form (1), then the 2016 pair products in
(2) assume at most

\[
 3^6-2^6=665
\tag{3}
\]

distinct values.

*Proof.* Match the 64 eigenvalue occurrences bijectively to the 64 bit vectors
\(\epsilon\in\{0,1\}^6\), even when two displayed values coincide. A pair of
distinct slots corresponds to \(\epsilon\ne\delta\), and its product is

\[
 a^2\prod_{i=1}^6 u_i^{\epsilon_i+\delta_i}.
\tag{4}
\]

The exponent vector \(f=\epsilon+\delta\) lies in \(\{0,1,2\}^6\). It has
at least one coordinate equal to 1: if every coordinate were 0 or 2, then
\(\epsilon_i=\delta_i\) for every \(i\), contrary to the two slots being
distinct. Conversely, there are exactly \(3^6-2^6\) vectors in
\(\{0,1,2\}^6\) having at least one 1.

Equation (4) maps every pair product into this set of at most 665 displayed
monomials. Equal or zero `u_i`, a zero `a`, accidental multiplicative
relations, or repeated eigenvalues can only merge displayed values; they can
never create a new exponent class. Thus degeneracies and zero modes only make
the bound stronger. []

### [LEMMA 2 — derivative gcd counts collisions in characteristic zero]

For a monic degree-\(N\) polynomial \(P\) over a characteristic-zero field,
let \(r\) be its number of distinct roots in an algebraic closure, with roots
counted only once for this number. Then

\[
 \deg\gcd(P,P')=N-r.
\tag{5}
\]

*Proof.* Write

\[
 P(z)=\prod_{\nu\in\mathcal D}(z-\nu)^{m_\nu},
\]

where \(\mathcal D\) is the set of distinct roots. At a root \(\nu\),

\[
 P'(z)=(z-\nu)^{m_\nu-1}
 \left(m_\nu h(z)+(z-\nu)h'(z)\right),\qquad h(\nu)\ne0.
\]

Characteristic zero makes \(m_\nu\ne0\), so the parenthesis is nonzero at
\(\nu\). Hence the common-root multiplicity is exactly \(m_\nu-1\), including
when \(\nu=0\). Therefore

\[
 \gcd(P,P')=\prod_{\nu\in\mathcal D}(z-\nu)^{m_\nu-1},
\]

and its degree is \(\sum_\nu(m_\nu-1)=N-r\). The gcd degree does not change
on extending the base field to an algebraic closure. []

Applying Lemmas 1--2 to (2) gives the necessary condition

\[
 \boxed{\deg\gcd(C_2,C_2')\ge 2016-665=1351}
\tag{6}
\]

for every full six-mode Gaussian spectrum. This proof used neither eigenvalue
ordering nor positivity.

## 3. Integral model removes the modular denominator gap

Directly reducing a rational `C_2` is not automatically enough: even if its
displayed coefficients are reducible at `p`, the primitive integral
representative of its monic rational gcd can have leading coefficient divisible
by `p`, so that representative can lose degree after reduction. Avoiding only
the denominators of the *entries* of `R` does not rule out this
factor-denominator exceptional-prime gap.

### [LEMMA 3 — monic integral reduction is degree-monotone at every prime]

Let \(R\in\mathbb Q^{64\times64}\), let \(D>0\) clear all of its entry
denominators, and set \(M=DR\in\mathbb Z^{64\times64}\). If
\(\rho_i=D\lambda_i\), define

\[
 C_M(z)=\prod_{i<j}(z-\rho_i\rho_j).
\tag{7}
\]

Then:

1. \(C_M\in\mathbb Z[z]\) is monic of degree 2016;
2. its characteristic-zero gcd degree is the same collision count as for
   \(C_2\), because all pair products are rescaled by the nonzero \(D^2\);
3. for **every** prime \(p\),
   \[
   \deg\gcd_{\mathbb F_p}(\overline {C_M},\overline {C_M'})
   \ \ge\ 
   \deg\gcd_{\mathbb Q}(C_M,C_M').
   \tag{8}
   \]

*Proof.* The induced action \(\bigwedge^2M\) is an integer matrix. After
triangularizing \(M\) over an algebraic closure, its diagonal on the exterior
square consists of \(\rho_i\rho_j\) for \(i<j\); this remains true with
algebraic multiplicity even if `M` is not diagonalizable. Thus
\(C_M=\operatorname{charpoly}(\bigwedge^2M)\), proving (1).

Let \(g\) be the monic gcd of \(C_M\) and \(C_M'\) over \(\mathbb Q[z]\).
Since \(C_M\) is monic integral, Gauss's lemma gives
\(g\in\mathbb Z[z]\), monic, and \(C_M/g\in\mathbb Z[z]\). Division of the
integral derivative \(C_M'\) by the monic integral divisor `g` also has an
integral quotient. Reducing these two factorizations modulo any prime leaves
\(\bar g\) as a common divisor. Its leading coefficient is still 1, so its
degree cannot drop. This proves (8). The nonzero scaling from `R` to `M`
preserves equality relations among pair products, proving (2). []

Thus finite-field gcd degrees are **upper bounds** on the desired
characteristic-zero degree. Reduction may add collisions and increase the
modular gcd degree, but cannot remove a characteristic-zero common factor.

## 4. Exact computation

**[COMPUTATION — rational input and integer scale].** The producer imports the
same `build_R` used by e38; it does not introduce another transfer convention.
For the open `2x3` bonds at \(t=1/3\), it finds

\[
 e^{2K}=\frac53,\qquad \epsilon=1,\qquad
 D=8{,}968{,}066{,}875=3^{15}5^4.
\tag{9}
\]

The script forms `D*R` using `Fraction` arithmetic and asserts every one of
its 4096 entries is integral. The factor \(5^4\) is necessary: the bond sum can
be \(b=-7\), so \((b-\epsilon)/2=-4\) and the diagonal factor contains
\((3/5)^4\). This is why `D` is computed as the actual entry-denominator LCM,
not inferred from a shortcut.

For each of the exact primes

\[
 p_1=1{,}000{,}003,\qquad p_2=2{,}000{,}003,
\]

primality is checked by trial division. Both are greater than 2016, so every
Newton divisor \(1,\ldots,2016\) and 2 is invertible. The artifact records the
literal inequalities

\[
 64(p-1)^2<2^{63},\qquad 2016(p-1)^2<2^{63},
\]

before using `int64` matrix products or modular Newton dot products; otherwise
the implementation would use Python integers.

For \(\bar M=M\bmod p\), the calculation is:

\[
\begin{aligned}
 s_k&=\operatorname{tr}(\bar M^k), &&1\le k\le64,\\
 \chi(z)&=z^{64}+a_1z^{63}+\cdots+a_{64},
 & k a_k&=-\sum_{i=1}^k a_{k-i}s_i,\\
 s_k&=-\sum_{i=1}^{64}a_i s_{k-i}, &&65\le k\le4032,\\
 q_m&=\sum_{i<j}(\rho_i\rho_j)^m
      =\frac{s_m^2-s_{2m}}2, &&1\le m\le2016.
\end{aligned}
\tag{10}
\]

Newton identities applied to the \(q_m\) reconstruct the monic degree-2016
pair polynomial. An explicit finite-field Euclidean algorithm then computes its
gcd with the formal derivative.

The characteristic-polynomial recurrence is checked independently against
binary matrix powering at

\[
 k\in\{65,66,67,100,128,256,512,1000,1024,2016,2017,3000,4032\}.
\]

For every full-size case and both primes, the producer additionally verifies
monicity and degrees (64 and 2016), evaluates \(\chi(\bar M)=0\) by matrix
Horner evaluation, and verifies that the returned gcd divides both polynomial
inputs with zero remainder.

### [COMPUTATION — results]

| exact case | integral scale | \(\deg\gcd\) at \(p_1\) | \(\deg\gcd\) at \(p_2\) |
|---|---:|---:|---:|
| `2x3` open 3D layer | \(3^{15}5^4\) | **385** | **385** |
| open chain `n=6` (2D strip control) | \(3^{14}5^3\) | 1351 | 1351 |
| synthetic full six-mode Gaussian | 1 | 1351 | 1351 |

For the target layer, Lemma 3 and the first row imply

\[
 \deg\gcd_{\mathbb Q}(C_M,C_M')\le385,
 \qquad
 r\ge2016-385=1631>665,
\tag{11}
\]

where \(r\) is the characteristic-zero number of distinct pair products. This
contradicts (3), proving Theorem P.

The full target coefficient-sequence digests (ascending coefficient order) are:

| prime | `charpoly_sha256` | `pair_poly_sha256` | `gcd_sha256` |
|---:|---|---|---|
| 1000003 | `677746d77fa33ca65b41e43f880ba5e01018869da93d4e5cf1224f0bdf105291` | `bd22fcfb805e33e7234f793aa404a0407c2a48dd90880c38e9cf3734d1d6ea9b` | `3b789be9bb92f9bcf869fe7568bb174081c1bf6235d365c74a3664f69cfad15d` |
| 2000003 | `dec7b3034342daf667bbd4fbaa29348d835d86ce8a4e5b0c5c74427d610333bc` | `5d430b9853de5f1bd9cb9bc90572cd530eac3cf3655b68ebae72d68cd7499fc2` | `1cfccf3dfb176eda99003121a7531c4312606204617f87aa341c3bc8c67b02b2` |

The JSON records the corresponding control digests, primes, all guard values,
all recurrence spot checks, degrees, timings, and all 48 producer checks.

### [COMPUTATION — independent controls]

1. The synthetic integer matrix is a unimodular conjugate of the known spectrum
   \(\{\prod_i u_i^{\epsilon_i}\}\) with
   \((u_1,\ldots,u_6)=(2,3,5,7,11,13)\). Unique factorization gives exactly
   665 distinct pair products over \(\mathbb Q\), hence exact gcd degree 1351.
   At both primes the producer reconstructs both its characteristic polynomial
   and its degree-2016 pair polynomial equal to the directly multiplied
   polynomials over the known eigenvalues, and obtains 1351.
2. The `n=6` open-chain 2D Ising strip control also gives 1351 at both primes.
   This is a consistency check, not an exact positive Gaussianity theorem from
   this calculation alone.
3. At both primes, two smaller independent ground truths (a deterministic random
   integer `6x6` matrix and the rationally built open chain `n=3`) reconstruct a
   pair polynomial equal to the determinant-interpolated
   \(\operatorname{charpoly}(\bigwedge^2M)\). Their direct exterior-square
   power traces also equal the pair-power formula in (10).

The standalone test independently rebuilds the target rational matrix rather
than importing either producer, repeats the degree-385 computation at \(p_1\),
and checks Newton signs, the pair-power formula, recurrence spots through 4032,
polynomial degrees/monicity, the chain control, and a separate exact four-mode
threshold control (`120 - (3^4-2^4)=55`). It then compares only recomputed
content fields/digests against the stored artifact.

## 5. Proof of Theorem P

**[THEOREM P, proof].** Suppose the `2x3` layer spectrum at \(t=1/3\) had the
full six-mode form (1). Lemmas 1--2 would force its characteristic-zero
pair-polynomial gcd degree to be at least 1351. By the exact integral-model
Lemma 3, the modular gcd degree at \(p_1=1000003\) would also be at least 1351.
The exact computation gives 385 instead. This contradiction rules out (1) for
`R`; similarity and nonzero scalar rescaling transfer the conclusion to the
physical layer operator. []

## 6. Why matching modular values are not overinterpreted

**[UNRESOLVED — exact characteristic-zero collision count].** Both selected
primes yield 385. This agreement is an independent implementation cross-check,
but it does **not** prove that the characteristic-zero gcd degree equals 385:
finite-field reduction can add common factors. The proved statement is only the
sufficient upper bound

\[
\deg\gcd_{\mathbb Q}(C_M,C_M')\le385
\]

and the resulting lower bound \(r\ge1631\). Exact equality is unnecessary for
the no-go theorem and is not claimed.

## 7. What this certificate does not decide

**[UNRESOLVED — parity and larger-space restrictions].** This pair-product
bound is for one full 64-value six-mode family. It does not itself decide:

* an even/odd parity-sector half, or a pair of separately parameterized
  parity-sector spectra;
* a restriction/submultiset of a Gaussian on a larger Hilbert space;
* a different coupling, an interval of couplings, a larger layer, or a
  thermodynamic-limit statement.

Those are different quantifiers. In particular, this certificate must not be
read as an all-\(K\), all-size, or parity-projected no-go. It is an exact,
basis-independent obstruction only to the stated full six-mode spectrum at the
single rational point \(t=1/3\).
