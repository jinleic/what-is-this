# Ordering-free pair-product obstruction for the open `2x5` layer

Artifacts: `experiments/e136_pair_product_2x5.py`,
`tests/test_pair_product_2x5.py`, and
`results/spectral/pair_product_2x5.json`.

**[THEOREM] Status.** The results below are an ordering-free finite-layer
certificate and an all-but-finitely-many-couplings corollary.  They concern the
full 1024-dimensional open `2x5` layer transfer spectrum, not a parity sector,
a larger Gaussian space restricted to the layer, or a thermodynamic limit.

## 1. Statements

**[THEOREM P10 — exact-coupling full-spectrum no-go].** At

\[
 t=\tanh(K^*/2)=\frac13,\qquad e^{2K}=\frac53,
\]

the open `2x5` layer transfer operator does not have a full ten-mode
subset-product spectrum

\[
 \left\{a\prod_{j=1}^{10}u_j^{\epsilon_j}:
   \epsilon\in\{0,1\}^{10}\right\}
\tag{1}
\]

as a multiset of 1024 eigenvalue slots, even if \(a,u_j\) are arbitrary complex
numbers.

**[THEOREM G10 — generic full-spectrum no-go].** There is a finite set
\(E\subset\mathbb C\), with

\[
 |E|\le 51{,}576{,}065{,}587{,}200,
\tag{2}
\]

such that every positive real \(t\notin E\) has an open `2x5` layer transfer
spectrum that is not of the form (1).  The bound is deliberately crude, and
this theorem neither locates nor identifies any member of \(E\).

**[COMPUTATION] Decisive finite-field values.** For the actual-LCM integral
model \(M=DR\) defined in §3, the exact modular calculations give

\[
 \deg\gcd\bigl(\operatorname{charpoly}(\bigwedge^2M),
                  \partial_z\operatorname{charpoly}(\bigwedge^2M)\bigr)
 =124{,}921
\tag{3}
\]

over both \(\mathbb F_{1{,}000{,}003}\) and
\(\mathbb F_{2{,}000{,}003}\).  Since \(124{,}921<465{,}751\), this is the
certificate behind Theorem P10.

## 2. Pair-product collision lemmas

Let \(\lambda_1,\ldots,\lambda_{2^{10}}\) denote eigenvalue *slots*, repeated
by algebraic multiplicity, and define

\[
 C_2(z)=\prod_{1\le i<j\le 2^{10}}(z-\lambda_i\lambda_j).
\tag{4}
\]

**[LEMMA 1 — Gaussian pair-product count].** If the slot multiset has the full
\(n\)-mode form
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

**[LEMMA 3 proof].** This is exactly the `e107_pair_product_obstruction.py`
actual-LCM monic integral convention \(R_{\rm int}=DR\), here applied at
dimension 1024.  The exterior-square matrix is integral and has the products
\(D^2\lambda_i\lambda_j\), \(i<j\), as its eigenvalue slots.  Its monic
characteristic-zero gcd is monic integral by Gauss's lemma; the two
corresponding monic quotient polynomials are integral as well.  Reducing those
factorizations modulo \(p\) leaves that monic gcd as a common divisor, whose
degree cannot drop.  Nonzero common rescaling preserves pair-product equality
multiplicities. []

**[LEMMA 4 — the ten-site Gaussian floor].** For \(n=10\),

\[
 \binom{1024}{2}=523{,}776,
 \quad 3^{10}-2^{10}=58{,}025,
 \quad 523{,}776-58{,}025=465{,}751.
\tag{9}
\]

Thus a full ten-mode spectrum of the form (1) requires
\(\deg\gcd(C_2,C_2')\ge465{,}751\).

**[LEMMA 4 proof].** Substitute \(n=10\) into Lemmas 1–2. []

## 3. Canonical rational matrix and its actual LCM

**[COMPUTATION] Canonical input.** The computation uses the existing
`e38_gaussianity_certificate.build_R` convention,

\[
 R=P_t\operatorname{diag}(q^{e_s})P_t,
 \qquad (P_t)_{xy}=t^{d_H(x,y)},
 \qquad q=\frac{1+t^2}{2t}.
\tag{10}
\]

At \(t=1/3\), \(q=5/3\).  The open `2x5` graph has 13 bonds (8 horizontal,
5 vertical), degree profile \(2^4 3^6\), is bipartite, and its canonical bond
sums give \(\epsilon=1\) and \(e_s\in[-7,6]\).

**[LEMMA 5 — exact split-integral construction].** Put

\[
 \widetilde P=3^{10}P_t,
 \qquad \widetilde d_s=5^{e_s+7}3^{6-e_s},
 \qquad W=\widetilde P\operatorname{diag}(\widetilde d_s)\widetilde P.
\tag{11}
\]

Then \(W\in\mathbb Z^{1024\times1024}\) and

\[
 R=\frac{W}{N_0},
 \qquad
 N_0=3^{20}\cdot5^7\cdot3^6=3^{26}5^7
   =198{,}583{,}267{,}838{,}203{,}125.
\tag{12}
\]

If \(G=\gcd(\gcd_{ij}W_{ij},N_0)\), the actual LCM of all reduced entry
denominators of \(R\) is \(D=N_0/G\), and \(M=DR=W/G\) is integral.

**[LEMMA 5 proof].** Since \(3^{10}t^h=3^{10-h}\) and
\((5/3)^e\,5^7 3^6=5^{e+7}3^{6-e}\), multiplying (10) by the displayed common
denominator gives (11)–(12); both exponents \(e+7\) and \(6-e\) are
nonnegative exactly because \(e\in[-7,6]\).  The reduced denominator of
\(W_{ij}/N_0\) is \(N_0/\gcd(W_{ij},N_0)\); the LCM of these divisors of
\(N_0\) is \(N_0/\gcd(\gcd_{ij}W_{ij},N_0)\).  Finally
\(D(W/N_0)=W/G\). []

**[COMPUTATION] Actual LCM result.** The exact content scan finds \(G=1\), so

\[
 D=N_0=3^{26}5^7=198{,}583{,}267{,}838{,}203{,}125.
\tag{13}
\]

The producer checks direct `Fraction` row sums on a first-row stride and 16
deterministically selected off-row entries.  The standalone verifier instead
reconstructs every one of the \(1024\times1024\) numerators \(W_{ij}\) by an
exact pure-Python integer route (each entry an individual
\(\texttt{math.sumprod}\) of Python integers), recomputes \(G\) by a full scan
with no early break, derives the actual LCM a second way as the running LCM of
every entry's reduced denominator \(N_0/\gcd(W_{ij},N_0)\), anchors 23
deterministic entries to true `Fraction` arithmetic, and independently
obtains (13).

**[COMPUTATION] Exact implementation guard.** To avoid a dense `Fraction`
triple loop in the producer, the Hamming kernel is built by a vectorized
popcount, and
\(L=\widetilde P\operatorname{diag}(\widetilde d)=2^{26}L_{\rm hi}+L_{\rm lo}\)
with \(0\le L_{\rm lo}<2^{26}\), so
\(W=2^{26}(L_{\rm hi}\widetilde P)+L_{\rm lo}\widetilde P\).  The code asserts
before both int64 products that every product's 1024-term integer sum is
below \(2^{63}\).  Therefore this is exact integer arithmetic, not a
numerical approximation.

## 4. Exact finite-field reconstruction

**[COMPUTATION] Prime and arithmetic conditions.** Both declared primes are
prime by trial division, are odd, avoid \(3\) and \(5\) (so 2 and \(G\)
invert), and exceed the polynomial degree

\[
 N=\binom{1024}{2}=523{,}776<p_1=1{,}000{,}003,
\tag{14}
\]

so every Newton divisor \(k\le N\) and every derivative exponent \(1\le k\le
N\) is a unit mod \(p\).  The trace-extension index \(2N=1{,}047{,}552\)
exceeds \(p_1\), but it is an *array index*, never a residue; this is the one
place the wave-13 requirement \(p>2N\) is loosened to the exact requirement
\(p>N\).  The code records and asserts

\[
 523{,}776\,(p-1)^2<2^{63},
 \qquad 1024\,(p-1)^2<2^{53},
 \qquad p+(523{,}777)(p-1)^2<2^{63}
\tag{15}
\]

for each prime.  The first bound makes every modular Newton dot product exact
in int64, the second makes every binary64 BLAS matrix-product partial sum an
exactly representable integer, and the third bounds the lazy-division drift.

**[LEMMA 6 — binary64 trace path is exact].** The producer's float64 matrix
products compute the same modular products as exact integer multiplication for
this input.

**[LEMMA 6 proof].** Every matrix factor is reduced to \([0,p)\) before a
product.  Every product is an integer below \((p-1)^2\), and every nonnegative
partial sum is below the second bound in (15), hence below \(2^{53}\).  Such
integers and their sums are represented exactly in IEEE binary64; changing the
summation order cannot exceed the same nonnegative total.  After each product
the code asserts that the BLAS output round-trips through int64 unchanged, and
one full int64 matrix product plus int64 binary matrix powers audit the whole
route. []

**[COMPUTATION] Reconstruction chain.** For each prime, the script computes

\[
\begin{aligned}
 s_k&=\operatorname{tr}(M^k), &&1\le k\le1024,\\
 \chi_M(z)&=z^{1024}+a_1z^{1023}+\cdots+a_{1024},
 & k a_k&=-\sum_{i=1}^k a_{k-i}s_i,\\
 s_k&=-\sum_{i=1}^{1024}a_i s_{k-i}, &&1025\le k\le1{,}047{,}552,\\
 q_m&=\frac{s_m^2-s_{2m}}2, &&1\le m\le523{,}776.
\end{aligned}
\tag{16}
\]

A second Newton reconstruction applied to \(q_m\) produces the monic
\(C_2=\operatorname{charpoly}(\bigwedge^2M)\) of degree 523,776.  Both Newton
loops and the trace recurrence dot contiguous int64 slices against
fill-from-the-end buffers, so no hidden strided copies occur; the code checks
\(\chi_M(M)=0\) by float64 Horner, recurrence values against int64 and float64
binary powers at clamped spot indices (including both endpoints \(N\) and
\(2N\)), and divisibility of both \(C_2\) and \(C_2'\) by the returned monic
gcd.

**[COMPUTATION] Lazy Euclidean exactness.** In the producer's polynomial
division, each window update subtracts \(\ell\,Q(z)\), where \(Q\) is the
current divisor, in int64 before one final tail reduction modulo \(p\).  One
entry can be touched by up to \(\deg Q+1\) successive windows, so its drift is
bounded by \(p+(\deg Q+1)(p-1)^2\le p+523{,}777\,(p-1)^2<2^{63}\), as
asserted.  Reading the current leading coefficient modulo \(p\) performs the
identical division in \(\mathbb F_p[z]\).  The standalone verifier uses strict
modular reduction after every update, a pure-Python ascending-coefficient
reference division on random inputs, and gets the same polynomials, hashes,
and gcd degrees.

## 5. Finite results and controls

**[COMPUTATION] Full exact result table.** The table gives finite-field gcd
degrees, not an asserted equality to the characteristic-zero degree.

| case | \(D\) | pair slots | Gaussian maximum | gcd at \(p_1\) | gcd at \(p_2\) | distinct-pair lower bound |
|---|---:|---:|---:|---:|---:|---:|
| open `2x3`, \(t=1/3\) | \(3^{15}5^4\) | 2,016 | 665 | 385 | 385 | 1,631 |
| open `3x3`, \(t=1/3\) | \(3^{24}5^6\) | 130,816 | 19,171 | 59,641 | 59,641 | 71,175 |
| open `2x5`, \(t=1/3\) | \(3^{26}5^7\) | 523,776 | 58,025 | **124,921** | **124,921** | **398,855** |
| open chain \(n=10\), \(t=1/3\) | \(3^{24}5^5\) | 523,776 | 58,025 | 465,751 | 465,751 | 58,025 |

**[COMPUTATION] Regression controls.** The `2x3` and `3x3` rows replay the
wave-12 and wave-13 gcd values at both primes, and their characteristic
polynomial, pair-polynomial, and gcd SHA-256 coefficient digests match the
stored artifacts exactly.  The open-chain \(n=10\) control — the same
523,776-slot trace/Newton/gcd path as the decisive case — reaches *exactly*
the ten-mode Gaussian floor \(465{,}751\) at both primes, with Euclid step
count \(58{,}025=3^{10}-2^{10}\); this checks the full pipeline at the
decisive size without asserting a new all-coupling positive theorem.

**[THEOREM P10 proof].** Lemma 3 and the `2x5` row of the table give
\(\deg\gcd_{\mathbb Q}(C_2,C_2')\le124{,}921\).  Lemma 2 therefore gives at
least

\[
 523{,}776-124{,}921=398{,}855
\tag{17}
\]

distinct pair products in characteristic zero.  This exceeds the Gaussian
maximum \(58{,}025\) in Lemma 4, contradicting (1).  The canonical \(R\)
differs from the physical transfer operator only by similarity and a nonzero
global scalar, which preserve collision multiplicities; this proves the stated
transfer-layer claim. []

## 6. Generic-\(t\) specialization/resultant argument

**[LEMMA 7 — polynomial representative and clearing exponent].** For the open
`2x5` layer, \(e\in[-7,6]\), and the clearing exponent must satisfy
\(c\ge\max(-e_{\min},e_{\max})=7\).  For \(c=7\),

\[
 M(t)=(2t)^{14}q(t)^{7}R(t)\in\mathbb Z[t]^{1024\times1024}.
\tag{18}
\]

Every entry of \(M(t)\) has \(t\)-degree at most 47.

**[LEMMA 7 proof].** Inside (10), the diagonal term after applying (18) is

\[
 (2t)^{14}q^{7+e}=(2t)^{7-e}(1+t^2)^{7+e},
\tag{19}
\]

which is an integer polynomial for every \(-7\le e\le6\) — note that both
\(c-e\ge0\) and \(c+e\ge0\) are needed, so the wave-13 choice \(c=F=6\)
fails at \(e=-7\) where \(F+e=-1\); the asymmetric exponent range of `2x5`
forces \(c=E=7\).  The two kernel entries contribute at most
\(d_H(i,s)+d_H(s,j)\le20\).  The degree in (19) is
\((7-e)+2(7+e)=21+e\le27\), giving \(20+27=47\). []

Define

\[
 A(t,z)=\operatorname{charpoly}(\bigwedge^2M(t))
       =z^{523776}+a_1(t)z^{523775}+\cdots+a_{523776}(t).
\tag{20}
\]

**[LEMMA 8 — coefficient and factor bounds].** Every exterior-square entry
has degree at most 94, and \(\deg_t a_k\le94k\).  The same indexed bound
holds for monic factors of \(A\) in \(\mathbb Q[t][z]\) and for the quotient
of \(\partial_zA\) by such a factor.

**[LEMMA 8 proof].** An exterior-square entry is a difference of two products
of degree-at-most-47 entries.  The coefficient \(a_k\) is a sum of
\(k\)-by-\(k\) minors.  For the factor statement, substitute \(z=t^{94}y\),
divide by the appropriate leading power of \(t\), and use monic factorization
in the UFD \(\mathbb Q[t^{-1}]\), exactly as in the wave-12 proof.  Monic
long division propagates the indexed degree bound to the derivative quotient.
[]

**[LEMMA 9 — specialization bound].** Let

\[
 G(t,z)=\gcd_{\mathbb Q(t)[z]}(A,\partial_zA),\qquad g=\deg_zG.
\tag{21}
\]

Then \(g\le124{,}921\).

**[LEMMA 9 proof].** At \(t=1/3\), (18) is

\[
 M(1/3)=\frac{2^{14}}{3^{47}}\,(DR),
\tag{22}
\]

where \(D\) is (13).  The nonzero scalar preserves the characteristic-zero
pair-collision gcd degree, which is at most 124,921 by §5.  Because \(A\) and
its monic generic gcd are monic in \(z\), both the gcd and its complementary
factor lie in \(\mathbb Q[t][z]\), so specialisation at \(t=1/3\) leaves a
monic common divisor of the same degree.  Hence \(g\) cannot exceed the
specialized degree. []

**[LEMMA 10 — finite exceptional set].** Put \(A_0=A/G\),
\(B_0=(\partial_zA)/G\), and

\[
 \rho(t)=\operatorname{Res}_z(A_0,B_0).
\tag{23}
\]

Then \(\rho\ne0\), and any \(t_0\) with \(\rho(t_0)\ne0\) satisfies
\(\gcd(A(t_0,z),\partial_zA(t_0,z))=G(t_0,z)\).

**[LEMMA 10 proof].** The quotients are coprime over \(\mathbb Q(t)\), so
their resultant is nonzero.  Resultants commute with specialization; when the
specialized quotients remain coprime, reintroducing the common monic factor
proves the equality. []

**[LEMMA 11 — effective exceptional bound].** The nonzero \(\rho\) in (23)
has

\[
 \deg_t\rho\le2\cdot94\cdot523{,}776\cdot523{,}775
 =51{,}576{,}065{,}587{,}200.
\tag{24}
\]

**[LEMMA 11 proof].** If \(n_0=\deg_zA_0\), the Sylvester determinant has
bidegree \((n_0-1,n_0)\) in the coefficients of \((A_0,B_0)\).  Lemma 8
bounds the indexed coefficient degrees by \(94k\), hence its degree is at
most \(2\cdot94\,n_0(n_0-1)\), and \(n_0\le523{,}776\) gives (24). []

**[THEOREM G10 proof].** Let \(E\) be the roots of \(\rho\).  Lemmas 9–11
show that every \(t>0\) outside \(E\) has gcd degree at most 124,921, hence
at least 398,855 distinct pair products.  Lemma 4 excludes (1), and (24)
bounds \(|E|\). []

## 7. Resource accounting and scope

**[COMPUTATION] Predeclared budgets.** Every compute gate in the producer
reads `time.process_time()` (CPU seconds), never wall clock and never
`signal.alarm`; wall-clock seconds are recorded as reference only.  The
predeclared per-stage CPU budgets are 900 s for construction, 900 s for
traces, 600 s for characteristic reconstruction, 900 s for the trace
recurrence, 1,200 s for pair-Newton reconstruction, and 4,500 s for the
Euclidean gcd, with a 9,000 s per-case CPU budget and a 6,000,000,000-byte
RSS cap.  A budget overrun aborts only its own case, recording the blocking
stage and its measured numbers; no stage was blocked in the recorded run.

**[COMPUTATION] Observed resource use.** The decisive `2x5` case measured,
per prime \((p_1,p_2)\): traces 13.580/13.304 CPU s, characteristic Newton
plus Horner 12.095/11.945 CPU s, trace recurrence plus audited binary powers
60.934/58.998 CPU s, pair Newton 44.241/43.480 CPU s, and Euclidean gcd
414.357/412.458 CPU s, with wall totals of 627.0/587.8 s per pipeline and
4.29 CPU s for the integral build.  The \(n=10\) chain control ran the same
sized pipeline in 275.3/276.3 wall seconds per prime with gcd stage
115.371/114.390 CPU s.  The full producer completed in 1,870.051 wall seconds
(1,661.913 CPU seconds) with measured peak RSS 278,069,248 bytes, and its
33-item check ledger passed in full.  These are observed measurements; the
preceding budgets are preflight declarations, not observed limits.

**[COMPUTATION] Standalone-verifier budgets.** The independent verifier gates
on `time.process_time()` with a 10,800 s overall CPU budget and a separate
5,400 s local CPU clock for each strict-gcd stage (probed by a test that a
freshly exhausted local clock fires), so prior case construction cannot
consume a later gcd budget.  It rebuilds the decisive `2x5` case and the
chain \(n=10\) control from the raw rational construction at \(p_1\), and the
`2x3`/`3x3` control rows from full `Fraction`-object matrices, then compares
every stored digest row by row.

**[UNRESOLVED] Deliberate exclusions.** The calculation does not assert that
124,921 is the characteristic-zero gcd degree; it is only the rigorous
modular upper bound needed here.  It does not enumerate the roots of \(\rho\),
decide any exceptional positive coupling, address parity-projected Gaussian
sectors, or imply a solution of the infinite-volume three-dimensional Ising
model.  Layers with six or more columns (dim \(\ge4096\),
\(\binom{4096}{2}\approx8.4\) million pair slots) remain open on this
pipeline.
