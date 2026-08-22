# Ordering-free pair-product obstruction for the open `3x3` layer

Artifacts: `experiments/e122_pair_product_3x3.py`,
`tests/test_pair_product_3x3.py`, and
`results/spectral/pair_product_3x3.json`.

**[THEOREM] Status.** The results below are an ordering-free finite-layer
certificate and an all-but-finitely-many-couplings corollary.  They concern the
full 512-dimensional open `3x3` layer transfer spectrum, not a parity sector,
a larger Gaussian space restricted to the layer, or a thermodynamic limit.

## 1. Statements

**[THEOREM P9 — exact-coupling full-spectrum no-go].** At

\[
t=\tanh(K^*/2)=\frac13,\qquad e^{2K}=\frac53,
\]

the open `3x3` layer transfer operator does not have a full nine-mode
subset-product spectrum

\[
\left\{a\prod_{j=1}^{9}u_j^{\epsilon_j}:
  \epsilon\in\{0,1\}^{9}\right\}
\tag{1}
\]

as a multiset of 512 eigenvalue slots, even if \(a,u_j\) are arbitrary complex
numbers.

**[THEOREM G9 — generic full-spectrum no-go].** There is a finite set
\(E\subset\mathbb C\), with

\[
 |E|\le 2{,}874{,}932{,}766{,}720,
\tag{2}
\]

such that every positive real \(t\notin E\) has an open `3x3` layer transfer
spectrum that is not of the form (1).  The bound is deliberately crude, and
this theorem neither locates nor identifies any member of \(E\).

**[COMPUTATION] Decisive finite-field values.** For the actual-LCM integral
model \(M=DR\) defined in §3, the exact modular calculations give

\[
 \deg\gcd\bigl(\operatorname{charpoly}(\bigwedge^2M),
                  \partial_z\operatorname{charpoly}(\bigwedge^2M)\bigr)
 =59{,}641
\tag{3}
\]

over both \(\mathbb F_{1{,}000{,}003}\) and
\(\mathbb F_{2{,}000{,}003}\).  Since \(59{,}641<111{,}645\), this is the
certificate behind Theorem P9.

## 2. Pair-product collision lemmas

Let \(\lambda_1,\ldots,\lambda_{2^n}\) denote eigenvalue *slots*, repeated
by algebraic multiplicity, and define

\[
 C_2(z)=\prod_{1\le i<j\le 2^n}(z-\lambda_i\lambda_j).
\tag{4}
\]

**[LEMMA 1 — Gaussian pair-product count].** If the slot multiset has the
full \(n\)-mode form
\(\{a\prod_j u_j^{\epsilon_j}:\epsilon\in\{0,1\}^n\}\), then its
\(\binom{2^n}{2}\) distinct-slot pair products assume at most

\[
 3^n-2^n
\tag{5}
\]

distinct values.

**[LEMMA 1 proof].** A pair of distinct bit strings \(\epsilon\ne\delta\)
has coordinatewise exponent sum
\(\epsilon+\delta\in\{0,1,2\}^n\), with at least one coordinate equal to
one.  There are exactly \(3^n-2^n\) such exponent vectors.  Coincident
parameters, zero parameters, repeated eigenvalues, and accidental
multiplicative relations can only merge the displayed monomials, so the bound
allows arbitrary complex parameters. []

**[LEMMA 2 — derivative gcd counts collisions].** For a monic degree-\(N\)
polynomial \(P\) over characteristic zero with \(r\) distinct roots,

\[
 \deg\gcd(P,P')=N-r.
\tag{6}
\]

**[LEMMA 2 proof].** Writing
\(P=\prod_{\nu\in\mathcal D}(z-\nu)^{m_\nu}\), the common multiplicity at
\(\nu\) of \(P\) and \(P'\) is \(m_\nu-1\) in characteristic zero.
Summing gives \(\sum_\nu(m_\nu-1)=N-r\). []

**[LEMMA 3 — monic integral reduction].** If \(R\in\mathbb Q^{d\times d}\),
\(D\) clears every entry denominator, and \(M=DR\in\mathbb Z^{d\times d}\),
then

\[
 A_M(z)=\operatorname{charpoly}(\bigwedge^2M)\in\mathbb Z[z]
\tag{7}
\]

is monic and, for every prime \(p\),

\[
 \deg\gcd_{\mathbb Q}(A_M,A_M')
 \le
 \deg\gcd_{\mathbb F_p}(\overline A_M,\overline A_M').
\tag{8}
\]

**[LEMMA 3 proof].** The exterior-square matrix is integral and has the
products \(D^2\lambda_i\lambda_j\), \(i<j\), as its eigenvalue slots.  Its
monic characteristic-zero gcd is monic integral by Gauss's lemma; the two
corresponding monic quotient polynomials are integral as well.  Reducing those
factorizations modulo \(p\) leaves that monic gcd as a common divisor, whose
degree cannot drop.  Nonzero common rescaling preserves pair-product equality
multiplicities. []

**[LEMMA 4 — the nine-site Gaussian floor].** For \(n=9\),

\[
 \binom{512}{2}=130{,}816,
 \quad 3^9-2^9=19{,}171,
 \quad 130{,}816-19{,}171=111{,}645.
\tag{9}
\]

Thus a full nine-mode spectrum of the form (1) requires
\(\deg\gcd(C_2,C_2')\ge111{,}645\).

**[LEMMA 4 proof].** Substitute \(n=9\) into Lemmas 1--2. []

## 3. Canonical rational matrix and its actual LCM

**[COMPUTATION] Canonical input.** The computation uses the existing
`e38_gaussianity_certificate.build_R` convention,

\[
 R=P_t\operatorname{diag}(q^{e_s})P_t,
 \qquad (P_t)_{xy}=t^{d_H(x,y)},
 \qquad q=\frac{1+t^2}{2t}.
\tag{10}
\]

At \(t=1/3\), \(q=5/3\).  The open `3x3` graph has 12 bonds, is bipartite,
and its canonical bond sums give \(\epsilon=0\) and
\(e_s\in[-6,6]\).

**[LEMMA 5 — exact split-integral construction].** Put

\[
 \widetilde P=3^9P_t,
 \qquad \widetilde d_s=5^{e_s+6}3^{6-e_s},
 \qquad W=\widetilde P\operatorname{diag}(\widetilde d_s)\widetilde P.
\tag{11}
\]

Then \(W\in\mathbb Z^{512\times512}\) and

\[
 R=\frac{W}{N_0},
 \qquad N_0=3^{24}5^6=4{,}412{,}961{,}507{,}515{,}625.
\tag{12}
\]

If \(G=\gcd(\gcd_{ij}W_{ij},N_0)\), the actual LCM of all reduced entry
denominators of \(R\) is \(D=N_0/G\), and \(M=DR=W/G\) is integral.

**[LEMMA 5 proof].** Since
\(3^9t^h=3^{9-h}\) and
\((5/3)^e3^65^6=5^{e+6}3^{6-e}\), multiplying (10) by the displayed common
denominator gives (11)--(12).  The reduced denominator of \(W_{ij}/N_0\) is
\(N_0/\gcd(W_{ij},N_0)\); the LCM of these divisors of \(N_0\) is
\(N_0/\gcd(\gcd_{ij}W_{ij},N_0)\).  Finally \(D(W/N_0)=W/G\). []

**[COMPUTATION] Actual LCM result.** The exact content scan finds \(G=1\),
so

\[
 D=N_0=3^{24}5^6=4{,}412{,}961{,}507{,}515{,}625.
\tag{13}
\]

The producer checks direct `Fraction` row sums on a first-row stride and 16
deterministically selected off-row entries.  The standalone verifier instead
constructs every one of the 512-by-512 rational entries, computes the LCM from
all of their denominators, and independently obtains (13).

**[COMPUTATION] Exact implementation guard.** To avoid a dense `Fraction`
triple loop in the producer, write
\(L=\widetilde P\operatorname{diag}(\widetilde d)=2^{26}L_{\rm hi}+L_{\rm lo}\)
with \(0\le L_{\rm lo}<2^{26}\), and compute
\(W=2^{26}(L_{\rm hi}\widetilde P)+L_{\rm lo}\widetilde P\).  The code
asserts before both int64 products that every product's 512-term integer sum
is below \(2^{63}\).  Therefore this is exact integer arithmetic, not a
numerical approximation.

## 4. Exact finite-field reconstruction

**[COMPUTATION] Prime and arithmetic conditions.** Both declared primes are
prime by trial division, exceed \(2\cdot130{,}816\), and are coprime to
\(3\cdot5\).  The code records and asserts

\[
 130{,}816(p-1)^2<2^{63},
 \qquad 512(p-1)^2<2^{53}
\tag{14}
\]

for each prime.  The first bound makes every modular Newton dot product and
lazy-reduction drift exact in int64.  The second makes every binary64 BLAS
matrix-product partial sum an exactly representable integer.

**[LEMMA 6 — binary64 trace path is exact].** The producer's float64 matrix
products compute the same modular products as exact integer multiplication for
this input.

**[LEMMA 6 proof].** Every matrix factor is reduced to \([0,p)\) before a
product.  Every product is an integer below \((p-1)^2\), and every nonnegative
partial sum is below the second bound in (14), hence below \(2^{53}\).  Such
integers and their sums are represented exactly in IEEE binary64; changing the
summation order cannot exceed the same nonnegative total.  After each product,
the integer result is reduced modulo \(p\).  The implementation additionally
compares the first full binary64 product to an int64 product and checks several
traces using all-int64 binary matrix powers. []

**[COMPUTATION] Reconstruction chain.** For each prime, the script computes

\[
\begin{aligned}
 s_k&=\operatorname{tr}(M^k), &&1\le k\le512,\\
 \chi_M(z)&=z^{512}+a_1z^{511}+\cdots+a_{512},
 & k a_k&=-\sum_{i=1}^k a_{k-i}s_i,\\
 s_k&=-\sum_{i=1}^{512}a_i s_{k-i}, &&513\le k\le261{,}632,\\
 q_m&=\frac{s_m^2-s_{2m}}2, &&1\le m\le130{,}816.
\end{aligned}
\tag{15}
\]

A second Newton reconstruction applied to \(q_m\) produces the monic
\(C_2=\operatorname{charpoly}(\bigwedge^2M)\) of degree 130,816.  The code
checks \(\chi_M(M)=0\), recurrence values against int64 binary powers through
\(k=261{,}632\), and divisibility of both \(C_2\) and \(C_2'\) by the returned
monic gcd.

**[COMPUTATION] Lazy Euclidean exactness.** In the producer's polynomial
division, each window update subtracts \(\ell\,Q(z)\), where \(Q\) is the
current divisor, in int64 before one final tail reduction modulo \(p\).  At most 130,816 such updates touch an
entry, so its magnitude is bounded by
\(p+130{,}816(p-1)^2<2^{63}\), as asserted.  Reading the current leading
coefficient modulo \(p\) performs the identical division in \(\mathbb F_p[z]\).
The standalone verifier uses strict modular reduction after every update and
gets the same polynomials, hashes, and gcd degrees.

## 5. Finite results and controls

**[COMPUTATION] Full exact result table.** The table gives finite-field gcd
degrees, not an asserted equality to the characteristic-zero degree.

| case | \(D\) | pair slots | Gaussian maximum | gcd at \(p_1\) | gcd at \(p_2\) | distinct-pair lower bound |
|---|---:|---:|---:|---:|---:|---:|
| open `2x3`, \(t=1/3\) | \(3^{15}5^4\) | 2,016 | 665 | 385 | 385 | 1,631 |
| open `2x4`, \(t=1/3\) | \(3^{21}5^5\) | 32,640 | 6,305 | 9,329 | 9,329 | 23,311 |
| open `3x3`, \(t=1/3\) | \(3^{24}5^6\) | 130,816 | 19,171 | **59,641** | **59,641** | **71,175** |
| open chain `n=9`, \(t=1/3\) | \(3^{22}5^4\) | 130,816 | 19,171 | 111,645 | 111,645 | 19,171 |

**[COMPUTATION] Regression controls.** The `2x3` and `2x4` gcd values replay
the wave-12 values at both primes.  Their characteristic-polynomial,
pair-polynomial, and gcd SHA-256 coefficient digests also match the wave-12
artifact exactly.  The open-chain `n=9` control reaches exactly the
nine-mode Gaussian floor at both primes; this checks the same 130,816-slot
trace/Newton/gcd path without asserting a new all-coupling positive theorem.

**[THEOREM P9 proof].** Lemma 3 and the `3x3` row of the table give
\(\deg\gcd_{\mathbb Q}(C_2,C_2')\le59{,}641\).  Lemma 2 therefore gives at
least

\[
 130{,}816-59{,}641=71{,}175
\tag{16}
\]

distinct pair products in characteristic zero.  This exceeds the Gaussian
maximum 19,171 in Lemma 4, contradicting (1).  The canonical `R` differs from
the physical transfer operator only by similarity and a nonzero global scalar,
which preserve collision multiplicities; this proves the stated transfer-layer
claim. []

## 6. Generic-\(t\) specialization/resultant argument

**[LEMMA 7 — polynomial representative].** For the open `3x3` layer,
\(e\in[-6,6]\), and

\[
 M(t)=(2t)^{12}q(t)^6R(t)\in\mathbb Z[t]^{512\times512}.
\tag{17}
\]

Every entry of \(M(t)\) has \(t\)-degree at most 42.

**[LEMMA 7 proof].** Inside (10), the diagonal term after applying (17) is

\[
 (2t)^{12}q^{6+e}=(2t)^{6-e}(1+t^2)^{6+e},
\tag{18}
\]

which is an integer polynomial for every \(-6\le e\le6\).  The two kernel
entries contribute at most \(d_H(i,s)+d_H(s,j)\le18\).  The degree in (18) is
\((6-e)+2(6+e)=18+e\le24\), giving \(18+24=42\). []

Define

\[
 A(t,z)=\operatorname{charpoly}(\bigwedge^2M(t))
       =z^{130816}+a_1(t)z^{130815}+\cdots+a_{130816}(t).
\tag{19}
\]

**[LEMMA 8 — coefficient and factor bounds].** Every exterior-square entry
has degree at most 84, and \(\deg_t a_k\le84k\).  The same indexed bound
holds for monic factors of \(A\) in \(\mathbb Q[t][z]\) and for the quotient
of \(\partial_zA\) by such a factor.

**[LEMMA 8 proof].** An exterior-square entry is a difference of two products
of degree-at-most-42 entries.  The coefficient \(a_k\) is a sum of
\(k\)-by-\(k\) minors.  For the factor statement, substitute \(z=t^{84}y\),
divide by the appropriate leading power of \(t\), and use monic factorization
in the UFD \(\mathbb Q[t^{-1}]\), exactly as in the wave-12 proof.  Monic
long division propagates the indexed degree bound to the derivative quotient.
[]

**[LEMMA 9 — specialization bound].** Let

\[
 G(t,z)=\gcd_{\mathbb Q(t)[z]}(A,\partial_zA),\qquad g=\deg_zG.
\tag{20}
\]

Then \(g\le59{,}641\).

**[LEMMA 9 proof].** At \(t=1/3\), (17) is

\[
 M(1/3)=\frac{2^{12}}{3^{42}}\,(DR),
\tag{21}
\]

where \(D\) is (13).  The nonzero scalar preserves the characteristic-zero
pair-collision gcd degree, which is at most 59,641 by §5.  Because \(A\) and
its monic generic gcd are monic in \(z\), both the gcd and its complementary
factor lie in \(\mathbb Q[t][z]\), so specialisation at \(t=1/3\) leaves a
monic common divisor of the same degree.  Hence \(g\) cannot exceed the
specialized degree. []

**[LEMMA 10 — finite exceptional set].** Put \(A_0=A/G\),
\(B_0=(\partial_zA)/G\), and

\[
 \rho(t)=\operatorname{Res}_z(A_0,B_0).
\tag{22}
\]

Then \(\rho\ne0\), and any \(t_0\) with \(\rho(t_0)\ne0\) satisfies
\(\gcd(A(t_0,z),\partial_zA(t_0,z))=G(t_0,z)\).

**[LEMMA 10 proof].** The quotients are coprime over \(\mathbb Q(t)\), so
their resultant is nonzero.  Resultants commute with specialization; when the
specialized quotients remain coprime, reintroducing the common monic factor
proves the equality. []

**[LEMMA 11 — effective exceptional bound].** The nonzero \(\rho\) in (22)
has

\[
 \deg_t\rho\le2\cdot84\cdot130{,}816\cdot130{,}815
 =2{,}874{,}932{,}766{,}720.
\tag{23}
\]

**[LEMMA 11 proof].** If \(n_0=\deg_zA_0\), the Sylvester determinant has
bidegree \((n_0-1,n_0)\) in the coefficients of \((A_0,B_0)\).  Lemma 8
bounds the indexed coefficient degrees by \(84k\), hence its degree is at
most \(2\cdot84n_0(n_0-1)\), and \(n_0\le130{,}816\) gives (23). []

**[THEOREM G9 proof].** Let \(E\) be the roots of \(\rho\).  Lemmas 9--11
show that every \(t>0\) outside \(E\) has gcd degree at most 59,641, hence at
least 71,175 distinct pair products.  Lemma 4 excludes (1), and (23) bounds
\(|E|\). []

## 7. Resource accounting and scope

**[COMPUTATION] Predeclared walls.** The producer declared stage walls of 300 s
for construction, 240 s for traces, 120 s for characteristic reconstruction,
240 s for recurrence, 600 s for pair Newton reconstruction, and 1,200 s for
the Euclidean gcd, with a 7,500,000,000-byte RSS cap.

**[COMPUTATION] Observed resource use.** The final artifact records a
per-stage RSS snapshot as well as duration.  The completed `3x3` layer took
57.012 s at \(p_1\) and 57.451 s at \(p_2\).  Its measured stage times were,
respectively, 1.103/1.111 s for traces, 0.868/0.884 s for
characteristic-Newton plus Horner, 23.950/24.398 s for recurrence plus int64
spot powers, 5.673/5.714 s for pair Newton, and 25.418/25.345 s for the gcd
stage, including 5.818/5.792 s of exact divisibility checks.  The stage-RSS
snapshots never exceeded 140,034,048 bytes for these two target-prime
pipelines.  The full producer completed in 215.464 s with measured peak RSS
154,451,968 bytes.  These are observed measurements; the preceding walls are
preflight budgets, not observed limits.

**[COMPUTATION] Standalone-verifier walls.** The independent raw-Fraction
verifier has a predeclared 3,600 s total wall and a separate 900 s clock for
each strict-gcd stage (including its paired exact divisibility audits), so
prior case construction cannot consume a later gcd budget.  Its final run
completed in 529.63 s with peak RSS 139,902,976 bytes; its two decisive
strict-gcd pipelines took 107.113 s and 106.453 s.  These are observed
measurements, whereas the two walls are preflight budgets.

**[UNRESOLVED] Deliberate exclusions.** The calculation does not assert that
59,641 is the characteristic-zero gcd degree; it is only the rigorous modular
upper bound needed here.  It does not enumerate the roots of \(\rho\), decide
any exceptional positive coupling, address parity-projected Gaussian sectors,
or imply a solution of the infinite-volume three-dimensional Ising model.
