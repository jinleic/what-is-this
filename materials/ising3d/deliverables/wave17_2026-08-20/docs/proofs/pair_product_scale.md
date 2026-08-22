# Generic and scaled pair-product obstruction

Artifacts: `experiments/e118_pair_product_scale.py`,
`tests/test_pair_product_scale.py`, and
`results/spectral/pair_product_scale.json`.

**Status convention.** Every mathematical statement below is tagged
**[THEOREM]**, **[LEMMA]**, **[COMPUTATION]**, or **[UNRESOLVED]**.  Modular
computations are exact finite-field facts; they are used only in the
one-sided direction justified in §2.

## 1. Headline result

**[THEOREM] Generic full-spectrum no-go for the open `2x3` layer.** There
is a finite set \(E\subset\mathbb C\), with

\[
 |E|\le 422{,}472{,}960,
\]

such that, for every real \(t>0\) outside \(E\), the 64 eigenvalues of the
open `2x3` transfer layer do **not** have the full six-mode subset-product
form

\[
 \left\{a\prod_{j=1}^{6}u_j^{\epsilon_j}:
       \epsilon\in\{0,1\}^{6}\right\}.
\]

The theorem allows arbitrary complex \(a,u_j\); in particular it excludes
the positive fermionic-Gaussian subclass.  It is not an interval statement:
the exceptional values are finite but neither located nor evaluated here.

The proof is algebraic over \(\mathbb Q(t)\), not an extrapolation from a
grid.  The grid in §5 is only an independent deterministic computational
cross-check.

## 2. Pair products and the monic integral model

Let \(\lambda_1,\ldots,\lambda_{2^n}\) be eigenvalue *slots*, repeated by
algebraic multiplicity, and put

\[
 C_2(z)=\prod_{i<j}(z-\lambda_i\lambda_j).
\]

**[LEMMA] Gaussian pair-product bound (Lemma 1).** A full \(n\)-mode
subset-product spectrum has at most \(3^n-2^n\) distinct pair products.
Consequently it requires

\[
 \deg\gcd(C_2,C_2')\ge
 \binom{2^n}{2}-(3^n-2^n).
\]

*Proof.* A pair of distinct bit strings has coordinatewise exponent sum in
\(\{0,1,2\}^n\) and at least one coordinate equal to 1.  There are
\(3^n-2^n\) such vectors.  If a monic polynomial of degree \(N\) has \(r\)
distinct roots in characteristic zero, then
\(\deg\gcd(P,P')=N-r\), including repeated zero roots.  Apply this to
\(C_2\). []

**[LEMMA] Monic integral reduction (Lemma 2).** If \(R\in\mathbb Q^{d\times d}\)
and \(D\) is the actual LCM of every entry denominator, then
\(M=DR\in\mathbb Z^{d\times d}\) and

\[
 A_M(z)=\operatorname{charpoly}(\bigwedge^2M)\in\mathbb Z[z]
\]

is monic.  For every prime \(p\),

\[
 \deg\gcd_{\mathbb Q}(A_M,A_M')
 \le \deg\gcd_{\mathbb F_p}(\overline A_M,\overline A_M').
\]

*Proof.* The exterior-square matrix is integral and has the pair products of
the eigenvalue slots as its eigenvalues.  The monic characteristic-zero gcd,
and both corresponding quotients, are integral by Gauss's lemma.  Reducing
these monic factorizations modulo \(p\) cannot lower the divisor degree. []

**[COMPUTATION] Canonical construction.** Every finite case imports
`build_R` from `e38_gaussianity_certificate.py`, then forms the actual LCM
\(D\) from its `Fraction` entries.  It reconstructs \(A_M\) by the
wave-11 chain

\[
 \operatorname{tr}(M^k)\;\longrightarrow\;\chi_M
 \;\longrightarrow\;\operatorname{tr}(M^k),\ 1\le k\le2\binom d2
 \;\longrightarrow\;\operatorname{tr}((\bigwedge^2M)^m)
 \;\longrightarrow\; A_M
 \;\longrightarrow\;\gcd(A_M,A_M').
\]

Newton identities, the characteristic recurrence, and the Euclidean gcd are
all exact in each finite field.  Matrix Horner evaluation checks
Cayley--Hamilton, binary powers check recurrence samples through the full
required range, and the returned gcd is checked to divide both inputs.

## 3. Polynomial family for the generic theorem

For the open `2x3` layer write

\[
 q(t)=\frac{1+t^2}{2t},\qquad
 R(t)=P_t\operatorname{diag}(q(t)^{e_s})P_t.
\]

**[LEMMA] Polynomial representative (Lemma 3).** For this layer
\(e_s\in\{-4,-3,\ldots,3\}\).  Therefore

\[
 M(t)=(2t)^7q(t)^4R(t)=(2t)^3(1+t^2)^4R(t)
 \in\mathbb Z[t]^{64\times64},
\]

and every entry of \(M(t)\) has \(t\)-degree at most 26.

*Proof.* The canonical entry of \(P_t\) is \(t^{d_H}\), and two such
factors contribute degree at most \(12\).  For \(-4\le e\le3\),

\[
 (2t)^7q(t)^{4+e}=(2t)^{3-e}(1+t^2)^{4+e}
\]

is an integer polynomial of degree \(11+e\le14\).  This gives the stated
matrix and entry-degree bounds.  For \(t>0\), the scalar multiplying
\(R(t)\) is nonzero, so pair-product collision multiplicities are unchanged.
[]

Define

\[
 A(t,z)=\operatorname{charpoly}(\bigwedge^2M(t))
       =z^{2016}+a_1(t)z^{2015}+\cdots+a_{2016}(t).
\]

**[LEMMA] Coefficient and factor degree bounds (Lemma 4).** The entries of
\(\bigwedge^2M(t)\) have degree at most 52, and

\[
 \deg_t a_k\le52k. \tag{1}
\]

The same indexed bound holds for every monic factor of \(A\) in
\(\mathbb Q[t][z]\), and for the quotient obtained by monic-dividing
\(\partial_zA\) by such a factor.

*Proof.* An exterior-square entry is a difference of two products of entries
of \(M\), proving the 52 bound.  A \(k\)-by-\(k\) minor gives (1).
For the factor claim, substitute \(z=t^{52}y\) and divide by the relevant
leading power of \(t\).  The resulting monic polynomial has coefficients in
the UFD \(\mathbb Q[t^{-1}]\).  A monic factor over its fraction field is
again in that UFD, giving the indexed bound after undoing the substitution.
Finally, monic long division of \(\partial_zA\) propagates the same bound
coefficient by coefficient. []

## 4. Specialization and genericity

Let

\[
 G(t,z)=\gcd_{\mathbb Q(t)[z]}(A,\partial_zA),\qquad g=\deg_zG,
\]

where the gcd is monic.

**[LEMMA] One specialization bounds the generic gcd (Lemma 5).** For every
\(t_0\in\mathbb Q\),

\[
 g\le\deg_z\gcd(A(t_0,z),\partial_zA(t_0,z)).
\]

*Proof.* Since \(A\) is monic in \(z\), both \(G\) and \(A/G\) lie in
\(\mathbb Q[t][z]\).  Monic division of \(\partial_zA\) by \(G\) also
has its quotient in \(\mathbb Q[t][z]\).  Thus \(G(t_0,z)\) is still monic
of degree \(g\) and divides both specialized polynomials. []

**[COMPUTATION] Decisive specialization.** At \(t=1/3\), the actual-LCM
model has \(D=3^{15}5^4=8{,}968{,}066{,}875\).  At both
\(p=1{,}000{,}003\) and \(p=2{,}000{,}003\), the exact modular gcd degree
is 385.  Lemma 2 gives the characteristic-zero bound at this specialization,
and Lemma 5 consequently gives

\[
 g\le385. \tag{2}
\]

The polynomial \(M(1/3)\) above is a nonzero rational scalar multiple of
this actual-LCM model, so its characteristic-zero gcd degree is the same.

Put \(A_0=A/G\), \(B_0=(\partial_zA)/G\), and

\[
 \rho(t)=\operatorname{Res}_z(A_0,B_0)\in\mathbb Q[t].
\]

**[LEMMA] Finite exceptional set (Lemma 6).** \(\rho\ne0\).  For every
\(t_0\) satisfying \(\rho(t_0)\ne0\),

\[
 \gcd(A(t_0,z),\partial_zA(t_0,z))=G(t_0,z),
\]

so the specialized gcd degree is exactly \(g\).

*Proof.* \(A_0\) and \(B_0\) are coprime over \(\mathbb Q(t)\), hence
their resultant is nonzero.  Resultants commute with specialization.  If its
specialized value is nonzero, the two specialized quotients are coprime;
reintroducing the common monic factor \(G(t_0,z)\) proves the equality. []

**[LEMMA] Effective exceptional-set bound (Lemma 7).** The nonzero polynomial
\(\rho\) has degree at most

\[
 2\cdot52\cdot2016\cdot2015=422{,}472{,}960. \tag{3}
\]

*Proof.* If \(n_0=\deg_zA_0\), then \(\deg_zB_0=n_0-1\).  Lemma 4 bounds
an indexed coefficient of \(A_0\) by \(52i\) and an indexed coefficient of
\(B_0\) by \(52j\).  The Sylvester determinant is homogeneous of degree
\(n_0-1\) in the first coefficients and degree \(n_0\) in the second;
therefore its \(t\)-degree is at most
\(2\cdot52n_0(n_0-1)\), which is bounded by (3). []

**[THEOREM] Proof of §1.** Outside the roots of \(\rho\), Lemmas 5--6
and (2) give \(\deg\gcd(A,\partial_zA)\le385\).  Hence there are at
least

\[
 2016-385=1631>3^6-2^6=665
\]

distinct pair products.  Lemma 1 contradicts a full six-mode
subset-product spectrum.  Lemma 7 bounds the number of exceptional complex
parameters, and restricting to \(t>0\) proves the stated result. []

## 5. Finite certificates and controls

**[COMPUTATION] `2x3` coupling grid.** The predeclared union of the
repository points and the regular grid \(t=k/20\), \(4\le k\le10\), is

\[
 \left\{\frac15,\frac14,\frac3{10},\frac13,\frac7{20},
        \frac25,\frac9{20},\frac12\right\}.
\]

For every listed `2x3` layer coupling, both primes give gcd degree 385.  Thus
at every one of these named points there are at least \(2016-385=1631\)
distinct pair products, exceeding the Gaussian maximum 665.  These are
separate exact point certificates; their common value is not used in place of
the generic proof.

**[COMPUTATION] Larger layer.** For the open `2x4` layer at \(t=1/3\),
\(q=5/3\), there are 10 in-layer bonds, and the actual scale is

\[
 D=3^{21}5^5=32{,}688{,}603{,}759{,}375.
\]

The decisive exact values are

| layer | dim | pair slots | Gaussian maximum | required gcd floor | gcd at \(p_1\) | gcd at \(p_2\) | distinct-pair lower bound |
|---|---:|---:|---:|---:|---:|---:|---:|
| open `2x4`, \(t=1/3\) | 256 | 32640 | 6305 | 26335 | 9329 | 9329 | 23311 |

Lemma 2 gives \(\deg\gcd_{\mathbb Q}\le9329\), hence
\(23311>6305\).  Therefore this named `2x4` finite layer is not a full
eight-mode Gaussian subset-product spectrum.

**[COMPUTATION] Controls.** The corresponding open chains return exactly
the Gaussian floors at both primes: degree 1351 for every `n=6` coupling in
the preceding grid and degree 26335 for the `n=8`, \(t=1/3\) control.  They
exercise the identical trace/Newton/gcd pipeline at the same site counts and
provide the required free-fermion consistency controls.  The result artifact
records all coefficient digests, recurrence samples, scales, timings, and
prime-wise values.

## 6. Bounded `3x3` attempt and scope

**[UNRESOLVED] `3x3` resource wall.** The open `3x3` target has dimension
512, \(\binom{512}{2}=130816\) pair slots, Gaussian maximum
\(3^9-2^9=19171\), and required gcd floor 111645.  The canonical exact
`build_R` route completed its rational/integral setup and entered direct
modular power traces at \(t=1/3\), but the predeclared 90-second wall fired
there after 90.058 seconds.  Recorded process peak RSS was 87,474,176 bytes,
well below the 7,500,000,000-byte cap.  No `3x3` modular gcd degree and no
`3x3` spectral conclusion are claimed.

**[UNRESOLVED] Exclusions.** The generic theorem does not identify or
exclude its finite exceptional set, does not decide parity-projected sectors
or restrictions of larger Gaussian spaces, and does not make a
thermodynamic-limit claim.  Matching two modular degrees is an implementation
cross-check only; it is never asserted to equal the characteristic-zero gcd
degree.
