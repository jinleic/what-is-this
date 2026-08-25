# A trace-seven resultant obstruction on the bipartite locus

Artifacts: `experiments/e232_trace_resultant.py`,
`results/spectral/trace_resultant.json`, and the independent verifier
`tests/test_trace_resultant.py`.

## 1. Outcome and scope

**[THEOREM — exact finite bipartite obstruction].** The physical open `2x3`
Ising layer transfer operator at

\[
t=\tanh(K^*/2)=\frac13,
\qquad e^{2K}=\frac{1+t^2}{2t}=\frac53,
\]

is not a full six-mode subset-product spectrum. The proof uses only the first
seven power traces. It does not order eigenvalues, form a pair polynomial,
compute eigenvalues, choose a Majorana basis, localize graph edges, or invoke a
varying-conjugator argument.

**[LEMMA — sharpness inside the named trace system].** After traces one through
three eliminate three sextic coefficients, the trace-four, trace-five, and
trace-six equations generate a proper ideal over `Q`. Adding trace seven makes
the exact ideal equal to `(1)`. Thus trace seven is the first decisive row of
this particular consecutive low-trace elimination.

**[COMPUTATION — controls].** A synthetic rational six-mode spectrum satisfies
all four remaining equations exactly. The physical open six-site chain retains
a proper ideal through trace seven, so the test does not reject the
free-fermion control.

**[SCOPE — extended below].** The pointwise certificate above is one finite graph at one
exact physical coupling. Sections 6--9 extend the same seven-trace necessity system over
every \(0<t<1\) and dispose its three algebraic branches. The resulting theorem is still
only for the finite open `2x3` layer: it is not an all-size bipartite-grid theorem or an
exact thermodynamic solution.

## 2. Full subset products force one mode polynomial

Let `W` be an invertible `64 x 64` matrix whose eigenvalue slots, including
multiplicity, have the full six-mode form

\[
 \lambda_S=a\prod_{i\in S}v_i,
 \qquad S\subseteq\{1,\ldots,6\},
 \qquad a,v_i\in\mathbb C^\times .                    \tag{2.1}
\]

Choose square roots `z_i^2=v_i` and put

\[
 g=a\prod_{i=1}^6 z_i,
 \qquad \mu_\varepsilon=\prod_{i=1}^6z_i^{\varepsilon_i},
 \qquad \varepsilon_i\in\{-1,+1\}.                    \tag{2.2}
\]

For the physical application `W=R` is positive definite. This removes a
potential root-of-unity ambiguity in the center: in any representation (2.1),
`a=lambda_emptyset>0` and `v_i=lambda_{\{i\}}/a>0`. Choose the positive
square roots `z_i`. Then `g>0`, while direct multiplication of all 64 slots
gives `g^64=det W`; hence `g` is the unique positive 64th root of the
determinant. It is therefore exactly the physical scalar center used in
Section 3, not that center times an unchecked root of unity.

Then the slots of `g^{-1}W` are exactly the `mu_epsilon`. Define centered power
traces

\[
 r_k=g^{-k}\operatorname{tr}(W^k).
\]

Put

\[
 \tau_i=z_i+z_i^{-1},
 \qquad u_i=\tau_i^2,
 \qquad Q(u)=\prod_{i=1}^6(u-u_i).                      \tag{2.3}
\]

**[LEMMA — Lucas factorization].** If `L_0(x)=2`, `L_1(x)=x`, and
`L_k(x)=xL_{k-1}(x)-L_{k-2}(x)`, then

\[
 r_k=\prod_{i=1}^6 L_k(\tau_i).                         \tag{2.4}
\]

**Proof.** `L_k(z+z^{-1})=z^k+z^{-k}` by the recurrence. Summing the slots
`mu_epsilon^k` independently over each sign gives the product in (2.4). `[]`

The needed Lucas polynomials factor as

\[
\begin{aligned}
L_1(\tau)&=\tau,\\
L_2(\tau)&=u-2,\\
L_3(\tau)&=\tau(u-3),\\
L_4(\tau)&=u^2-4u+2,\\
L_5(\tau)&=\tau(u^2-5u+5),\\
L_6(\tau)&=(u-2)(u^2-4u+1),\\
L_7(\tau)&=\tau(u^3-7u^2+14u-7).
\end{aligned}                                             \tag{2.5}
\]

Because the mode count is even, equations (2.3)--(2.5) imply the following
necessary identities:

\[
\begin{aligned}
 Q(0)&=r_1^2,\\
 Q(2)&=r_2,\\
 Q(3)&=r_3/r_1,\\
 \operatorname{Res}_u(u^2-4u+2,Q)&=r_4,\\
 \operatorname{Res}_u(u^2-5u+5,Q)&=r_5/r_1,\\
 \operatorname{Res}_u(u^2-4u+1,Q)&=r_6/r_2,\\
 \operatorname{Res}_u(u^3-7u^2+14u-7,Q)&=r_7/r_1.
\end{aligned}                                             \tag{2.6}
\]

The physical transfer representative below is positive definite, so `r_1` and
`r_2` are strictly positive and the displayed divisions are legitimate. Over a
general field the same constraints can instead be denominator-cleared.

**[LEMMA — resultant sign].** Every polynomial on the left of the last four
rows in (2.6) is monic. For monic `h`,
`Res(h,Q)=product_{h(alpha)=0}Q(alpha)`. Grouping first over the six roots
`u_i` gives precisely the product of `h(u_i)` in (2.5); the even degree six
removes the possible sign. `[]`

Thus seven traces of a putative subset-product spectrum must arise from one
monic sextic `Q`. Reciprocity alone imposes no such joint trace law; this is why
(2.6) survives the bipartite blind locus of the trace-versus-inverse invariant.

## 3. Exact rational physical representative

For the open `2x3` graph let `B=sum_(ij) Z_iZ_j`, with seven bonds, and let
`b_s` be its diagonal value in spin state `s`. At rational `t`, set

\[
 P_t(s,s')=t^{d_H(s,s')},
 \qquad q=\frac{1+t^2}{2t},
 \qquad R=P_t\,\operatorname{diag}
       \left(q^{(b_s-1)/2}\right)P_t.                    \tag{3.1}
\]

All `b_s` are odd, so `R` is rational. The symmetric physical sandwich `S`
satisfies

\[
 S=(1-t^2)^{-6}e^K R.                                    \tag{3.2}
\]

The matrix `P_t` is invertible for `0<t<1`, and the diagonal in (3.1) is
strictly positive, so `R` and `S` are positive definite. The latter is similar
to the usual layer transfer product. Therefore testing `R` is spectrally
equivalent up to one positive scalar. If `g=(det R)^{1/64}>0`, (3.2) and
`det S=1` give

\[
 z:=g^2=(1-t^2)^{12}q^{-1}.                              \tag{3.3}
\]

Consequently all inputs in (2.6) are rational combinations of raw traces
`p_k=tr(R^k)`:

\[
\begin{array}{c|c}
\text{quantity}&\text{exact expression}\\ \hline
r_1^2&p_1^2/z\\
r_2&p_2/z\\
r_3/r_1&p_3/(p_1z)\\
r_4&p_4/z^2\\
r_5/r_1&p_5/(p_1z^2)\\
r_6/r_2&p_6/(p_2z^2)\\
r_7/r_1&p_7/(p_1z^3).
\end{array}                                               \tag{3.4}
\]

At `t=1/3`, the producer clears (3.1) by the exact integer scale

\[
 D=8\,968\,066\,875
\]

and computes `tr((DR)^k)` for `k=1,...,7` by integer matrix multiplication.
The complete integers and the seven reduced rationals in (3.4) are stored in
`results/spectral/trace_resultant.json`; no floating-point value enters the
certificate.

## 4. Elimination certificate

Write

\[
 Q(u)=u^6+a_5u^5+a_4u^4+a_3u^3+a_2u^2+a_1u+a_0.          \tag{4.1}
\]

The first row of (2.6) fixes `a_0=r_1^2`. Let

\[
\begin{aligned}
 C_2&=r_2-a_0-64-32a_5-16a_4-8a_3,\\
 C_3&=r_3/r_1-a_0-729-243a_5-81a_4-27a_3.
\end{aligned}
\]

The `Q(2)` and `Q(3)` rows then give

\[
 a_2=C_3/3-C_2/2,
 \qquad a_1=3C_2/2-2C_3/3.                              \tag{4.2}
\]

Substituting (4.2) into the last four rows of (2.6) leaves four polynomials in
`Q[a_5,a_4,a_3]`, of total degrees

\[
 2,\ 2,\ 2,\ 3.                                         \tag{4.3}
\]

**[COMPUTATION — exact certificate].** Over `Q`, grevlex Groebner bases for the
first one, two, and three polynomials are proper, with basis sizes `1,3,6`.
For all four polynomials the reduced basis is exactly

\[
 \{1\}.                                                   \tag{4.4}
\]

Hence the four equations have no common point even over the algebraic closure.
Under assumption (2.1), the actual six values `u_i` would produce the sextic
(2.3) and therefore a common point, contradicting (4.4). This proves the
finite bipartite obstruction.

The independent verifier reconstructs (3.1) without importing the producer,
uses SymPy integer matrices instead of NumPy object matrices, computes each
cyclotomic norm as a determinant of multiplication in `Q[u]/(h)` instead of a
resultant, and obtains the same prefix verdict using lexicographic rather than
grevlex elimination.

## 5. Controls and method boundary

**[COMPUTATION — positive algebra control].** For

\[
 u_i\in\left\{\frac92,\frac{16}3,\frac{36}5,
                 \frac{64}7,\frac{144}{11},\frac{196}{13}\right\},
\]

constructing `Q=product_i(u-u_i)` makes all four residuals in (4.3) exactly
zero and recovers the stored `a_2,a_1`. This audits every Lucas/resultant sign.

**[COMPUTATION — physical free-fermion control].** Repeating the exact trace
construction for the open six-site chain at the same `t` leaves a proper ideal
through trace seven. This is a non-rejection control; the proper ideal alone is
not promoted to a sufficiency statement.

**[THEOREM — comparison with the prior pair-product obstruction].** The new
certificate uses seven `64 x 64` power traces and a three-variable elimination.
Its producer uses under `80 MB` peak RSS and under one process CPU second in the
recorded run. It is logically independent of the much larger unordered-pair
characteristic polynomial: neither theorem implies the other's certificate,
although both refute the same finite subset-product hypothesis at this point.

**[UNRESOLVED].** To turn this into an every-coupling bipartite theorem one
must control the physical-curve parameter dependence of the four-polynomial
ideal, for example by an explicit parameter resultant with a no-positive-root
certificate. One empty specialization proves the finite theorem and suggests
generic failure; by itself it does not remove all exceptional couplings.

## 6. Symbolic physical curve

Artifacts for the symbolic continuation are
`experiments/e238_trace_exceptional_set.py`,
`results/spectral/trace_exceptional_set.json`, and the independent verifier
`tests/test_trace_exceptional_set.py`.

**[THEOREM — finite exceptional set at every physical coupling].** For the
open `2x3` layer and every \(0<t<1\), a full six-mode subset-product spectrum
can occur only at one of three algebraic values of \(t\). There is exactly one
candidate in each of the three rational intervals displayed in Section 8.
The seven-trace system has an algebraic common point at every candidate, so
this theorem does **not** assert that the physical exceptional set is empty.

Keep

\[
 q=\frac{1+t^2}{2t}.
\]

The map \(t\mapsto q\) is a decreasing bijection from \((0,1)\) to
\((1,\infty)\). For a spin state \(s\), put
\(d_s=(b_s-1)/2\), and define the polynomial matrix

\[
 B(q)_{s,r}=q^{\,d_s+10-d_H(s,r)}.                       \tag{6.1}
\]

The exponent is nonnegative: on the open `2x3` graph,
\(-4\leq d_s\leq3\). Writing
\(H_k(q)=\operatorname{tr}(B(q)^k)\), cyclicity of trace and
\(P_t^2=(2t)^6(q^{6-d_H(s,r)})_{s,r}\) give

\[
 p_k=\operatorname{tr}(R^k)
     =(2t)^{6k}q^{-4k}H_k(q).                            \tag{6.2}
\]

Also,

\[
 z=(1-t^2)^{12}q^{-1}
   =(2t)^{12}(q^2-1)^6q^{-1}.                            \tag{6.3}
\]

Consequently the seven inputs to (2.6) are rational functions of \(q\)
alone. With \(\Delta=q^2-1\),

\[
\begin{aligned}
 r_1^2&=\frac{H_1^2}{q^7\Delta^6},&
 r_2&=\frac{H_2}{q^7\Delta^6},&
 \frac{r_3}{r_1}&=\frac{H_3}{H_1q^7\Delta^6},\\
 r_4&=\frac{H_4}{q^{14}\Delta^{12}},&
 \frac{r_5}{r_1}&=\frac{H_5}{H_1q^{14}\Delta^{12}},&
 \frac{r_6}{r_2}&=\frac{H_6}{H_2q^{14}\Delta^{12}},\\
 \frac{r_7}{r_1}&=\frac{H_7}{H_1q^{21}\Delta^{18}}.
\end{aligned}                                             \tag{6.4}
\]

Every \(H_k\) has positive integer coefficients. Thus all denominators used
in (6.4) are nonzero for \(q>1\). The producer obtains the \(H_k\) by
single-process multiplication over \(\mathbb Z[q]\); their degrees are

\[
13,\ 26,\ 39,\ 52,\ 65,\ 78,\ 91.                       \tag{6.5}
\]

This derives the physical-curve equations independently of any endpoint or
generic-finiteness statement.

## 7. A specialization-stable parameter norm

Let \(F_4,F_5,F_6,F_7\) be the four residuals obtained after (4.2).
Their total degrees in \((a_5,a_4,a_3)\) are \(2,2,2,3\). The homogeneous
quadratic parts of the first three residuals are independent of \(q\):

\[
\begin{aligned}
h_4={}&4a_3^2+56a_3a_4+280a_3a_5+188a_4^2
       +1816a_4a_5+4252a_5^2,\\
h_5={}&5a_3^2+75a_3a_4+390a_3a_5+275a_4^2
       +2800a_4a_5+6980a_5^2,\\
h_6={}&6a_3^2+84a_3a_4+432a_3a_5+276a_4^2
       +2700a_4a_5+6318a_5^2.
\end{aligned}                                             \tag{7.1}
\]

An exact grevlex basis of \((h_4,h_5,h_6)\) has leading monomials

\[
 a_3^4,\quad a_5a_3^2,\quad a_4a_3^2,\quad
 a_5^2,\quad a_5a_4,\quad a_4^2.                         \tag{7.2}
\]

Hence the homogeneous ideal has no projective zero and has Hilbert function
\(1,3,3,1,0,\ldots\). Its eight standard monomials are

\[
 {\cal B}=(1,a_3,a_3^2,a_3^3,a_5,a_5a_3,a_4,a_4a_3).
                                                                  \tag{7.3}
\]

**[LEMMA — stable rank-eight quotient].** At every \(q>1\),
\(\overline{\mathbb Q}[a_5,a_4,a_3]/(F_4,F_5,F_6)\) has dimension eight,
counting multiplicity, and (7.3) is a basis.

**Proof.** The three affine quadrics have the fixed leading forms (7.1).
Those forms have no common projective point by (7.2), so there are no roots
at infinity. Bézout therefore gives \(2^3=8\) affine roots counting
multiplicity. Lifting exact homogeneous-basis identities from the \(h_j\) to
the \(F_j\) reduces every nonstandard monomial to (7.3), without dividing by
a parameter polynomial. Thus (7.3) spans an eight-dimensional quotient and
is a basis. \(\square\)

Reduce \(bF_7\) to (7.3) for each \(b\in{\cal B}\), and let \(M_7(q)\) be
the resulting \(8\times8\) multiplication matrix. Define

\[
 E(q)=\operatorname{primnum}\det M_7(q).                 \tag{7.4}
\]

The exact cleared determinant has degree \(1307\). Cancelling its degree-336
gcd with the collected column denominator leaves a primitive, square-free,
dense polynomial

\[
 E\in\mathbb Z[q],\qquad \deg E=971.                     \tag{7.5}
\]

Its 972 coefficients are stored in ascending order in the JSON artifact,
with SHA-256
`5c7be98f19d7c75a53a1bb171496b00c85eaa1af6879073824698c41ad6b1c77`.

**[LEMMA — necessity and exact survival].** A common point of
\(F_4=F_5=F_6=F_7=0\) at \(q>1\) forces \(E(q)=0\). Conversely, every
\(q>1\) root of \(E\) gives an algebraic common point of all four equations.

**Proof.** By the stable rank-eight lemma, \(M_7\) is multiplication by
\(F_7\) in a finite-dimensional algebra. Its determinant vanishes exactly
when \(F_7\) is not a unit. In an Artinian algebra this is equivalent to
\(F_7\) lying in a maximal ideal, hence to a common algebraic point. The
denominators cancelled in (7.4) are nonzero for \(q>1\) by (6.4). \(\square\)

Pulling (7.5) back to the requested physical parameter gives the primitive
reciprocal polynomial

\[
 E_t(t)=\operatorname{prim}\left[
 (2t)^{971}E\!\left(\frac{1+t^2}{2t}\right)\right],
 \qquad \deg E_t=1942.                                  \tag{7.6}
\]

Its coefficient digest is
`9767345ab6ee92ad5269d18a296a9fad3b6ba0b60105e3ceba66a3d8d68bd889`.

## 8. Exact physical root isolation and disposition

**[COMPUTATION — exact interval certificate].** The subresultant gcd
\(\gcd(E,E')\) is \(1\). The coefficient sequence of \(E(1+x)\) has exactly
three sign variations, so Descartes' rule gives at most three roots for
\(q>1\), counted with multiplicity. Exact rational evaluation gives opposite
endpoint signs on each of the following three disjoint intervals:

\[
\begin{array}{c|c|c}
 &q\text{-interval}&(\operatorname{sgn}E(q_-),
                       \operatorname{sgn}E(q_+))\\ \hline
1&(1456799/10^6,\ 1821/1250)&(-,+)\\
2&(3198407/10^6,\ 399801/125000)&(+,-)\\
3&(2397237/250000,\ 9588949/10^6)&(-,+).
\end{array}                                               \tag{8.1}
\]

Thus there is exactly one root in each interval and no other root for
\(q>1\). Monotonicity of \(q=(1+t^2)/(2t)\), checked by exact rational
endpoint inequalities, maps them into

\[
\begin{aligned}
&397429/10^6<t<39743/10^5,\\
&160347/10^6<t<40087/250000,\\
&10457/200000<t<26143/500000.
\end{aligned}                                             \tag{8.2}
\]

The independent verifier runs in two bounded stages. Its elimination stage
reconstructs the norm by degree-filtered affine Macaulay reductions rather
than lifted H-basis division. Its root stage independently rebuilds the
integer binomial shift \(E(1+x)\), applies Descartes' rule, and checks all three
rational sign brackets and the exact \(q\)-to-\(t\) map. An earlier exact
SymPy Sturm/subresultant call was observed to exceed a 1,800-second wall after
the Macaulay stage had passed; it was replaced rather than reported as a
mathematical failure.

**[COMPUTATION — all-coupling chain control].** For the open six-site chain,
the same symbolic construction admits an explicit rational coefficient
branch \(Q(q)\). Exact substitution makes all four residuals
\(F_4,F_5,F_6,F_7\) identically zero for every \(q>1\). This is stronger than
the one-point non-rejection control in Section 5 and confirms that the
parameter norm does not reject the known free-fermion family.

**[COMPUTATION — bounded producer].** The recorded single-process run used
444.245132 process CPU seconds and 210124800 bytes peak RSS. Its preflight
wall was an \(8\times8\) polynomial determinant with cleared column-degree
maxima

\[
(102,156,186,240,156,186,156,186),
\]

whose determinant-degree bound is \(1368\). The producer refused a generic
four-variable lexicographic closure in advance because the fixed rank-eight
route gives a bounded exact certificate. The declared hard acceptance walls
were 900 process CPU seconds and \(2\) GiB RSS.

**[RESOLVED — superseded by e242 and Section 9].** Each interval in (8.2) contains a
genuine algebraic solution of the seven-trace coefficient system. Section 9 reconstructs
the unique coefficient sextic at each root and proves that none has six mode roots
\(u_i\geq4\). Thus the three values survive the polynomial trace equations but are not
physical full six-mode subset-product spectra.
The hash-pinned e238 artifact remains the record of the finite-exceptional-set stage; the
semialgebraic disposition is a separate e242 artifact.

## 9. Semialgebraic disposition of the three exceptional branches

The disposition artifacts are
`experiments/e242_trace_branch_disposition.py`,
`results/spectral/trace_branch_disposition.json`, and
`tests/test_trace_branch_disposition.py`. They consume the certified \(E(q)\) and its
three isolating intervals from e238; they do not recompute the degree-971 determinant.

### 9.1 A unique coefficient sextic at each exceptional root

Fix one of the three real roots \(\alpha>1\) of \(E\), with its real embedding selected
by the stored 1000-bit dyadic subinterval. Let

\[
 {\cal A}_\alpha=
 \overline{\mathbb Q}[a_5,a_4,a_3]/(F_4,F_5,F_6).
\]

The stable-basis lemma gives \(\dim {\cal A}_\alpha=8\), with basis (7.3). Let
\(M_\alpha\) be multiplication by \(F_7\). The producer clears each column by a
polynomial denominator and obtains

\[
 P_\alpha=M_\alpha D_\alpha .
\]

The cleared matrix has the same exact content digest as the e238 norm matrix. Exact
rational interval evaluation proves that all eight diagonal entries of \(D_\alpha\)
are nonzero at every branch, so \(P_\alpha\) and \(M_\alpha\) have the same image and
rank. Since \(E\) divides the cleared determinant, \(E(\alpha)=0\) gives
\(\det P_\alpha=0\).

Write

\[
 C_{j0}=(-1)^j\det(P\text{ with row }j\text{ and column }0\text{ deleted}).
\]

The exact cofactor degrees at \(j=0,1,4,6\) are respectively
\(1208,1211,1208,1209\). On all three dyadic intervals, a mean-value rational
enclosure excludes zero from \(C_{00}\). Hence

\[
 \operatorname{rank}P_\alpha=7,\qquad
 \dim {\cal A}_\alpha/(F_7)=1.                           \tag{9.1}
\]

A one-dimensional unital algebra over the coefficient field is that field itself, so
(9.1) is one reduced coefficient point. Row zero of \(\operatorname{adj}P_\alpha\)
annihilates the image. In the basis
\((1,a_3,a_3^2,a_3^3,a_5,a_5a_3,a_4,a_4a_3)\), normalization by its nonzero
constant coordinate gives

\[
 a_3=\frac{C_{10}}{C_{00}},\qquad
 a_5=\frac{C_{40}}{C_{00}},\qquad
 a_4=\frac{C_{60}}{C_{00}}.                             \tag{9.2}
\]

The remaining coefficients \(a_2,a_1,a_0\) are then reconstructed from (4.2) and
the exact trace targets, not inferred from numerical eigenvalues. Every division in
this chain has a rational interval denominator excluding zero.

### 9.2 Hermite signatures and the \(u=4\) boundary

For a monic sextic \(Q(u)\), let \(s_k\) be its Newton power sums and form the real
Hermite matrix

\[
 H(Q)_{ij}=s_{i+j},\qquad 0\leq i,j<6.                   \tag{9.3}
\]

All six leading principal minors are nonzero on every branch. Their signs, the
corresponding exact \(LDL^{T}\) pivot signs, and the Hermite signatures are

\[
\begin{array}{c|c|c|c}
\text{branch}&\operatorname{sgn}(\Delta_1,\ldots,\Delta_6)
 &\operatorname{sgn}(d_1,\ldots,d_6)&\operatorname{sig}H\\ \hline
1&(+,+,-,-,+,+)&(+,+,-,+,-,+)&2\\
2&(+,-,-,+,+,+)&(+,-,+,-,+,+)&2\\
3&(+,+,+,+,+,+)&(+,+,+,+,+,+)&6.
\end{array}                                               \tag{9.4}
\]

Hermite's theorem identifies \(\operatorname{sig}H\) with the number of distinct
real roots of the square-free sextic. Branches 1 and 2 therefore have only two real
roots and cannot be physical six-mode spectra.

Branch 3 has six distinct real roots, so its location relative to \(4\) matters.
Exact interval expansion gives the descending coefficient signs

\[
 \operatorname{sgn}Q(4+x)=(+,-,-,-,+,+,+).              \tag{9.5}
\]

The constant sign in (9.5) proves \(Q(4)\neq0\). The sequence has two sign
variations, so Descartes' rule bounds the number of positive roots of \(Q(4+x)\),
counted with multiplicity, by two. Thus at most two branch-3 roots satisfy \(u>4\);
six roots in \([4,\infty)\) are impossible.

### 9.3 Every-coupling finite-layer theorem

**[THEOREM — empty physical exceptional set on the open `2x3` layer].** For every
\(0<t<1\), the positive-definite open `2x3` Ising layer transfer spectrum is not a
full six-mode subset-product spectrum.

**Proof.** Section 8 confines any such spectrum to the three roots of \(E\).
Positivity of all 64 spectral slots forces every mode ratio \(v_i>0\). Therefore

\[
 u_i=\left(\sqrt{v_i}+\frac1{\sqrt{v_i}}\right)^2\geq4.  \tag{9.6}
\]

Equations (9.4)--(9.5) show that none of the three unique trace-system sextics has
six roots satisfying (9.6), a contradiction. \(\square\)

The producer reconstructs the quotient by a lifted H-basis. The independent verifier
uses degree-filtered affine Macaulay reduction, reproduces all eight clearing
denominators, the full cleared-matrix digest, and four cofactor polynomials, then
independently redoes the dyadic nesting, interval arithmetic, Hermite signs, and
Descartes bound. The producer used 814.373944 process CPU seconds and
186515456 bytes peak RSS.

This theorem is finite and full-spectrum only. It does not exclude a parity-sector
factorization, a larger auxiliary construction, or exceptional points on other
bipartite grids, and it gives no thermodynamic free energy, critical coupling, or
critical exponent.
