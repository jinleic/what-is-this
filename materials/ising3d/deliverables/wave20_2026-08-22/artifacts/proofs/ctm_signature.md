# Finite corner-transfer signature: exact 2D control and a 3D definition obstruction

## 1. Scope and claim discipline

[UNRESOLVED] This note does **not** solve the standard three-dimensional Ising model, decide an infinite-volume corner spectrum, prove non-integrability, or give an all-size CTM theorem.  It defines one finite graph kernel exactly, proves a non-vacuous two-dimensional product-spectrum control, evaluates one declared three-dimensional boundary regularization, and proves that the honest three-dimensional tensor does not itself determine a matrix spectrum.

[COMPUTATION] The only coupling used for the point certificates is the exact rational low-temperature edge parameter $q=e^{-2K}=1/2$; this local symbol $q$ is the repository variable $x$.  The benchmark $K_c=0.221654626$ is not used to choose, fit, or assess any operator.

## 2. The finite boundary kernel is defined before spectral testing

[LEMMA] For $d\geq2$ and $r\geq1$, let
\[
 Q_d(r)=\{0,1,\ldots,r\}^d
\]
with open nearest-neighbour bonds, and let $F_i=\{x_i=0\}$ be its $i$th coordinate face.  A bit equal to $1$ denotes spin $+1$, as in the repository convention.  For face assignments $\alpha_i\in\{\pm1\}^{F_i}$, define
\[
 C_{d,r}(q)[\alpha_1,\ldots,\alpha_d]
 =\begin{cases}
 \displaystyle\sum_{\substack{\sigma:Q_d(r)\to\{\pm1\}\\
                  \sigma|_{F_i}=\alpha_i}}
 q^{b(\sigma)},&\text{the assignments agree on every face overlap},\\
 0,&\text{otherwise},
 \end{cases}
\]
where $b(\sigma)$ is the number of unequal-spin bonds.  Every graph bond occurs once; no half-bond convention or square root of a weight is introduced.

[LEMMA] For $d=2$, the two faces are boundary arms and $C_{2,r}$ is a symmetric matrix.  The shared face intersection is one corner spin.  Global spin flip makes its fixed-$+$ and fixed-$-$ blocks permutation-similar, so the fixed-$+$ block is a spectrum-independent representative of that duplication.

[LEMMA] For $d=3$, the three faces meet pairwise along edges and $C_{3,r}$ is a rank-three tensor.  Coordinate permutations permute its three legs, including the induced reordering of face-state bits.  A rank-three tensor is not an endomorphism, so it has no characteristic polynomial until an additional contraction is supplied.

[COMPUTATION] The radius choices are minimal for the stated tests.  The $d=2,r=2$ fixed-corner block has dimension $2^2=4$, the first size at which a two-mode product relation constrains the spectrum.  The $d=3,r=1$ octant is one cube with one spin outside the union of its three coordinate faces; its three face legs each have dimension $2^4=16$, and exactly $2^7=128$ compatible tensor entries are nonzero.

[COMPUTATION] As an independent convention check, summing both fixed-corner blocks of the $d=2$ kernel gives the scaled open-$3\times3$ partition value $70146$, exactly matching `ising.exact_enumeration.dos_bonds`.  Summing all entries of the $d=3$ tensor gives the scaled open-$2\times2\times2$ value $36450$, again exactly matching the repository enumerator.

## 3. Logarithm-free four-slot criterion

[LEMMA] Let $M$ be a real symmetric positive-definite $4\times4$ matrix with four simple eigenvalue slots $\lambda_0,\ldots,\lambda_3$.  The eigenvalues of $\wedge^2M$ are the six products $\lambda_i\lambda_j$ with $i<j$.  Thus the monic polynomial
\[
 P_2(z)=\chi_{\wedge^2M}(z)
\]
records all unordered slot-pair products without ordering eigenvalues or taking logarithms.

[LEMMA] If the four slots can be labelled as
\[
 c\{1,x,y,xy\},
\]
then $c\cdot(cxy)=(cx)(cy)$, so $P_2$ has a repeated root and $\gcd(P_2,P_2')$ has positive degree.  Conversely, when the four slots are nonzero and simple, equality between two pair products sharing an index would force two eigenvalues to be equal.  Therefore a repeated pair product must use two disjoint complementary pairs; relabelling those pairs gives $c\{1,x,y,xy\}$.  For a simple positive four-slot spectrum, a degree-one pair gcd is consequently an exact non-vacuous two-mode subset-product certificate.

[LEMMA] Any additive two-bit level formula for logarithms,
$\log\lambda_{n_1,n_2}=a+n_1b+n_2c$, exponentiates to the same complementary-product identity.  Failure of that identity therefore rejects this additive-in-exponents signature without evaluating a logarithm.

## 4. Exact 2D control

[THEOREM] Let $A_2(q)$ be the fixed-$+$ radius-two block of $C_{2,2}(q)$.  Exact arithmetic in $\mathbb Z[q]$ gives
\[
 \det A_2(q)=p(q)^2,
 \qquad p(q)=q^6(1-q^2)^6,
\]
and
\[
 \chi_{\wedge^2A_2(q)}(z)=(z-p(q))^2R(z,q)
\]
for a polynomial $R\in\mathbb Z[q,z]$.  Two exact synthetic divisions have zero remainder, while a third division has a nonzero polynomial remainder.  This is a symbolic relation, not a numerical eigenvalue coincidence.

[COMPUTATION] At $q=1/2$, multiplying every matrix entry by the common denominator $2^{12}=4096$ gives a positive-definite integer matrix.  Its characteristic polynomial is
\[
 z^4-12825z^3+14033412z^2-598363200z+2176782336,
\]
its characteristic gcd with its derivative is $1$, and its leading principal minors are
\[
 2160,\quad1213056,\quad18895680,\quad2176782336.
\]

[COMPUTATION] The pair-product polynomial of the scaled matrix is
\[
\begin{aligned}
 z^6&-14033412z^5+7671831257664z^4-654981671517659136z^3\\
 &+16699906766455659823104z^2
 -66495657533778638407729152z\\
 &+10314424798490535546171949056,
\end{aligned}
\]
and its exact gcd with its derivative is $z-46656$.  Hence there are exactly five distinct values among the six pair-product slots.

[THEOREM] The simple positive $d=2,r=2,q=1/2$ spectrum is exactly relabellable as $c\{1,x,y,xy\}$.  This is a non-vacuous finite two-dimensional integrable control for the proposed Baxter-style multiplicative/additive signature, proved without logarithms.

## 5. Honest 3D tensor and one declared spectralization

[LEMMA] The honest radius-one three-dimensional object is $C_{3,1}(q)\in V^{\otimes3}$ with $\dim V=16$.  It is not silently flattened.  For a transparent finite comparison only, declare
\[
 h_{\mathrm{free}}(\gamma)=1
\]
on the third face, contract that leg, and fix the two-site intersection edge of the remaining faces to spin $+1$.  The two unfixed spins on each remaining face index a $4\times4$ matrix $B_{\mathrm{free}}(q)$.

[COMPUTATION] At $q=1/2$, with the same common scale $4096$, the declared matrix is
\[
 B_{\mathrm{free}}=
 \begin{pmatrix}
 336&240&240&768\\
 240&264&192&960\\
 240&192&264&960\\
 768&960&960&5376
 \end{pmatrix}.
\]
Its leading principal minors are
\[
 336,\quad31104,\quad2737152,\quad1719926784,
\]
so it is positive definite, and
\[
 \chi_B(z)=z^4-6240z^3+2306880z^2-158008320z+1719926784
\]
is squarefree.

[COMPUTATION] Its exact pair-product polynomial is
\[
\begin{aligned}
 z^6&-2306880z^5+984251990016z^4-84001120934952960z^3\\
 &+1692841359833818988544z^2
 -6824092786556505717473280z\\
 &+5087798221017014024540258304,
\end{aligned}
\]
and $\gcd(P_2,P_2')=1$.  All six slot-pair products are distinct.

[THEOREM] The declared free-third-face, fixed-plus-edge radius-one 3D regularization is not a nonzero two-mode subset-product spectrum at $q=1/2$, and it admits no additive two-bit logarithmic level formula.  This is an exact finite boundary-condition counterexample only.

## 6. Definition-level obstruction

[THEOREM] The graph tensor, positivity, global spin flip, and square-face symmetry do not select a unique matrix spectrum.  Besides $h_{\mathrm{free}}$, define the strictly positive covector
\[
 h_{\mathrm{biased}}(\gamma)
 =1+\mathbf 1\{\text{all four spins of }\gamma\text{ are aligned}\}.
\]
Both covectors are invariant under global spin flip and every square-face symmetry, and they are not proportional.

[LEMMA] For a matrix $M$ with nonzero trace, the rational quantity
\[
 I_2(M)=\frac{\dim(M)\operatorname{tr}(M^2)}{\operatorname{tr}(M)^2}
\]
is unchanged by multiplying $M$ by a nonzero scalar.  Distinct $I_2$ values therefore prove that two spectra are not related merely by normalization.

[COMPUTATION] The full $16\times16$ contractions from $h_{\mathrm{free}}$ and $h_{\mathrm{biased}}$ both have exact rank $16$, but
\[
 I_2(M_{\mathrm{free}})=\frac{1308424}{294849},\qquad
 I_2(M_{\mathrm{biased}})=\frac{20362024}{4296645}.
\]
Their fixed-plus-edge sectors likewise have
\[
 \frac{5959}{1690}\neq\frac{94939}{26450}.
\]
Thus even normalization cannot identify the two spectra.

[COMPUTATION] The spin-flip-invariant aligned-face covector, supported only on the all-plus and all-minus face states, gives a full contraction of rank $8$ instead of $16$.  Pair-product collisions created by such zero modes are therefore a boundary-rank artifact, not evidence of integrability.

[THEOREM] A matrix spectrum can be obtained only after adding a third-face boundary covector or a different contraction rule.  Symmetry and positivity leave spectrally inequivalent choices.  Consequently the exact free-face failure above cannot be promoted to an intrinsic CTM spectrum of $C_{3,1}$, and selecting a contraction because it has a favorable collision pattern would be circular.

## 7. Counterexamples and residual risk

[COMPUTATION] The squarefree $B_{\mathrm{free}}$ pair polynomial is a counterexample to the naive claim that the smallest free-face 3D corner automatically inherits the finite 2D two-mode product signature.

[THEOREM] The free and biased positive covectors are a counterexample to the claim that the stated finite graph and its symmetries determine a unique 3D corner matrix up to scalar normalization.

[COMPUTATION] The rank-eight aligned contraction is a counterexample to treating an abundance of pair collisions as non-vacuous without first excluding zero-mode multiplicity.

[UNRESOLVED] A future route could add and justify a specific contraction, prove compatibility across radii, and control its thermodynamic limit.  No such selection or limit theorem is supplied here.  The radius-one cube is finite and unusually small, so neither its failure nor any alternative finite success is an all-size or infinite-volume theorem.

[UNRESOLVED] The finite kernel used here deliberately counts every orthant bond once.  It is not asserted to equal a particular normalized infinite Baxter CTM, whose boundary splitting and limiting normalization would require an additional comparison theorem.

## 8. Reproduction

[COMPUTATION] The producer command is:

```bash
.venv/bin/python experiments/e208_ctm_certificate.py
```

[COMPUTATION] The standalone clean-room verifier command is:

```bash
.venv/bin/python tests/test_ctm_signature.py
```

[COMPUTATION] The verifier imports none of the producer modules.  It enumerates the 3D tensor by its $2^7$ compatible boundary assignments, constructs $\wedge^2M$ explicitly, and obtains characteristic polynomials by a Leibniz determinant rather than the producer's power-trace/Newton route.
