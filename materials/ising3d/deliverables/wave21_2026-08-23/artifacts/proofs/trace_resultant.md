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

**[UNRESOLVED].** This is one finite graph at one exact physical coupling. It
does not prove the every-coupling theorem for `2x3`, an all-size theorem for
bipartite grids, or an exact thermodynamic solution. In particular, a proper
finite trace ideal is only consistency with a subset-product spectrum, never a
sufficiency theorem.

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
