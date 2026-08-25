# Octahedral-star matchgate and the exact positivity barrier

Artifacts:

- producer: `experiments/e244_octahedral_matchgate.py`;
- clean-room verifier: `tests/test_octahedral_matchgate.py`;
- producer output: `results/integrability/octahedral_matchgate.json`.

## 1. Status and conventions

`[THEOREM]` and `[LEMMA]` below have self-contained proofs. `[COMPUTATION]`
means a bounded exact integer/`Fraction` enumeration. `[EXTERNAL]` is used only
for a cited standard identity or terminology. `[UNRESOLVED]` marks conclusions
that this result deliberately does not make.

Let the six boundary spins be \(\sigma_i\in\{-1,1\}\), let the eliminated
center spin be \(\tau\), and put \(v=\tanh K\), with \(0<v<1\). For a subset
\(A\subseteq[6]\), write

\[
 \chi_A(\sigma)=\prod_{i\in A}\sigma_i.
\]

There are two different Walsh objects in this note:

1. \(g(A)\) is the Walsh signature of the **weight** obtained by summing the
   center spin;
2. \(c_A\) is a Walsh coefficient of the **logarithm** of that weight.

The matchgate and tensor-rank statements concern \(g\). The alternating-sign
theorem in Section 4 concerns \(c_A\). Keeping these objects separate is
essential.

## 2. Exact star Walsh signature

**[LEMMA — star Walsh identity].** For every degree \(d\),

\[
\begin{aligned}
 W_d(\sigma)
 &=\sum_{\tau=\pm1}\exp\!\left(K\tau\sum_{i=1}^d\sigma_i\right)\\
 &=\cosh(K)^d\left\{\prod_{i=1}^d(1+v\sigma_i)
                  +\prod_{i=1}^d(1-v\sigma_i)\right\}\\
 &=2\cosh(K)^d\sum_{\substack{A\subseteq[d]\\ |A|\text{ even}}}
      v^{|A|}\chi_A(\sigma).
\end{aligned}
\]

**Proof.** Use \(e^{K\tau\sigma_i}=\cosh K(1+v\tau\sigma_i)\),
expand the product, and sum \(\tau^{|A|}\) over \(\tau=\pm1\). That sum is two
for even \(|A|\) and zero for odd \(|A|\). ∎

For \(d=6\), after omitting the common positive factor
\(2\cosh(K)^6\), the signature is therefore

\[
 g(A)=\begin{cases}
 v^{|A|},&|A|\text{ even},\\
 0,&|A|\text{ odd}.
 \end{cases}
 \tag{2.1}
\]

**[COMPUTATION].** At each of \(v=1/3\) and \(v=1/2\), the producer builds all
64 spin-basis values from the two center-spin products and performs the full
\(64\times64\) rational Walsh transform. Every entry agrees with (2.1). The
verifier reconstructs the center-spin sum independently and repeats the
transform.

## 3. Finite bipartite cycle-space factorization

Let \(G=(U\sqcup V,E)\) be a finite bipartite graph. For each \(u\in U\),
expanding its eliminated-center factor chooses an even subset
\(A_u\subseteq\delta(u)\). Summing a retained spin at \(w\in V\) kills a term
unless an even number of the chosen incident edges meet \(w\). Since every
edge has exactly one endpoint in \(U\), the tuple \((A_u)_{u\in U}\) is in
bijection with an edge set \(F\subseteq E\). Its two parity conditions say
exactly \(\partial F=\varnothing\). Thus:

**[LEMMA — cycle-space contraction].** Let
\(g_u(\sigma_{\delta(u)})=
\sum_{\substack{A_u\subseteq\delta(u)\\ |A_u|\text{ even}}}
v^{|A_u|}\chi_{A_u}(\sigma)\) be the normalized star signature. Then

\[
 2^{-|V|}\sum_{\sigma_V\in\{\pm1\}^{V}}
   \prod_{u\in U}g_u(\sigma_{\delta(u)})
 =\sum_{F\subseteq E:\ \partial F=\varnothing}v^{|F|}.
 \tag{3.1}
\]

Equivalently, the raw retained-spin sum has the factor \(2^{|V|}\) on
the right. Restoring the eliminated-star factors \(2\cosh(K)^{\deg u}\)
gives the standard high-temperature identity

\[
 Z_G(K)=2^{|U|+|V|}\cosh(K)^{|E|}
        \sum_{F:\partial F=\varnothing}v^{|F|}.
\]

This is a finite normal-factor-graph contraction: local factors impose the
parity constraints, edge variables select \(F\), and the explicit
\(2^{-|V|}\) averages the retained-spin sum.

**[COMPUTATION — controls].** The producer compares three exact routes: local
even-star contraction, direct edge-subset parity enumeration, and the
normalized spin sum. The verifier instead uses a boundary-state dynamic
program and a fresh spin sum. The controls are

| graph | edge-subset bound | cycle dimension | cycle polynomial |
|---|---:|---:|---:|
| \(C_4=K_{2,2}\) | \(2^4=16\) | 1 | \(1+v^4\) |
| \(K_{2,3}\) | \(2^6=64\) | 2 | \(1+3v^4\) |

The factorization is exact. It is not a claim that an arbitrary nonplanar
network contraction is one Pfaffian.

## 4. Alternating signs of every even log-Walsh sector

This section strengthens the degree-six sign calculation inherited from
`e234`. Let \(d\) be even, set \(S(\sigma)=\sum_{i=1}^d\sigma_i\), and define,
for nonempty \(A\subseteq[d]\),

\[
 c_A(K)=2^{-d}\sum_{\sigma\in\{\pm1\}^d}
            \chi_A(\sigma)\log\cosh(KS(\sigma)).
 \tag{4.1}
\]

Global spin flip kills every odd \(|A|\), so only even sectors remain.

**[EXTERNAL — standard product].** NIST DLMF 4.36.E2 gives

\[
 \cosh y=\prod_{n=0}^{\infty}
 \left(1+a_n y^2\right),\qquad
 a_n=\frac{4}{\pi^2(2n+1)^2}>0.
 \tag{4.2}
\]

Also, directly by differentiating in \(a y^2\) (and fixing the value at zero),

\[
 \log(1+a y^2)=\int_0^\infty
 e^{-s}\frac{1-e^{-a s y^2}}{s}\,ds.
 \tag{4.3}
\]

**[THEOREM — all even-star sectors].** Let \(d\) be even and let
\(\varnothing\ne A\subseteq[d]\) have \(|A|=2r\). For every real \(K\ne0\),

\[
 \operatorname{sign} c_A(K)=(-1)^{r+1}.
 \tag{4.4}
\]

Consequently, at every even degree \(d\ge6\), all two-spin log couplings are
strictly positive, all four-spin couplings are strictly negative, all
six-spin couplings are strictly positive, and the higher even sectors continue
to alternate.

**Proof.** Apply (4.3) to one factor of (4.2), with \(y=KS\). The constant
term in (4.3) drops because \(A\ne\varnothing\) and
\(\mathbb E\chi_A=0\). For \(b=a_nK^2s>0\), the Gaussian Fourier identity
gives

\[
 e^{-bS^2}=\frac1{\sqrt{4\pi b}}
 \int_{\mathbb R}e^{-t^2/(4b)}e^{itS}\,dt.
\]

Independence of the spins then yields

\[
\begin{aligned}
 \mathbb E\!\left[\chi_Ae^{itS}\right]
 &=\left(\mathbb E[\sigma e^{it\sigma}]\right)^{2r}
   \left(\mathbb E[e^{it\sigma}]\right)^{d-2r}\\
 &=(i\sin t)^{2r}(\cos t)^{d-2r}\\
 &=(-1)^r\sin(t)^{2r}\cos(t)^{d-2r}.
\end{aligned}
\]

Both exponents in the final integrand are even. The Gaussian-weighted integral
is therefore nonnegative and is strictly positive because its integrand is not
almost everywhere zero. Hence
\(\mathbb E[\chi_Ae^{-bS^2}]\) has strict sign \((-1)^r\).
Equation (4.3) contributes an additional minus sign, so every product factor
contributes strict sign \((-1)^{r+1}\). All factors have the same sign, proving
(4.4).

The interchanges are controlled explicitly. The Walsh expectation is finite.
The product-log sum converges uniformly on the finite spin state space because
\(\log(1+a_nK^2S^2)\le a_nK^2d^2\) and \(\sum_n a_n<\infty\).
For the outer integral near \(s=0\),
\(\mathbb E[\chi_A S^{2k}]=0\) for \(k<r\), since the polynomial degree
\(2k\) is smaller than \(|A|=2r\). Taylor's remainder therefore gives
\[
 \left|\mathbb E[\chi_Ae^{-a_nK^2sS^2}]\right|
 \le \frac{(a_nK^2sd^2)^r}{r!},\qquad 0<s\le1.
\]
After division by \(s\), the right-hand side is a constant times
\(s^{r-1}\), which is integrable at zero because \(r\ge1\).
For \(s\ge1\) the expectation has absolute value at most one, so the
\(e^{-s}ds/s\) integral is absolutely convergent. This also justifies carrying
the strictly signed Gaussian expectation through the integral and the product
sum. The assumption that \(d\) is even is essential to this proof:
for odd \(d\), the remaining power of \(\cos t\) is odd and has no fixed sign.
No odd-degree sign theorem is claimed.
∎

**[COMPUTATION — exact finite controls].** No logarithm is evaluated
numerically. For every nonempty even sector at \(d=2,4,6\), the producer
enumerates the integer multiplicities of
\(\log\cosh(K|d-2j|)\) before the factor \(2^{-d}\). The verifier reconstructs
those multiplicities independently as

\[
 \sum_{p=0}^{|A|}\sum_{q=0}^{d-|A|}
 (-1)^p\binom{|A|}{p}\binom{d-|A|}{q}
 [\,|d-2(p+q)|=m\,].
\]

The top-sector rows are

| \(d\) | \(2^d c_{[d]}\) | sign certificate |
|---:|---|---|
| 2 | \(2\log\cosh(2K)\) | positive for \(K\ne0\) |
| 4 | \(2\log\cosh(4K)-8\log\cosh(2K)\) | negative |
| 6 | \(2\log\cosh(6K)-12\log\cosh(4K)+30\log\cosh(2K)\) | positive |

For the last two rows, put \(x=\cosh(2K)^2>1\). Then

\[
 e^{8c_{[4]}}=\frac{2x-1}{x^2},\qquad
 x^2-(2x-1)=(x-1)^2>0,
\]

and

\[
 e^{32c_{[6]}}=\frac{x^8(4x-3)}{(2x-1)^6},
\]

with

\[
 x^8(4x-3)-(2x-1)^6=(x-1)^3Q_6(x),
\]

\[
 Q_6(1+y)=8+54y+144y^2+188y^3+120y^4+33y^5+4y^6>0.
\]

These finite polynomial certificates corroborate the analytic theorem; the
all-degree proof does not depend on them. No large \(d=8,10,12\) coefficient
payload is stored.

## 5. An explicit all-\(d\) Pfaffian lemma

For an ordered set of size \(d\), define the skew matrix

\[
 (J_d)_{ij}=\begin{cases}
 1,&i<j,\\
 -1,&i>j,\\
 0,&i=j.
 \end{cases}
\]

**[LEMMA — all-\(d\) Pfaffian certificate].** For every \(r\ge0\),

\[
 \operatorname{Pf}(J_{2r})=1.
 \tag{5.1}
\]

Consequently, for every \(d\) and every even subset \(A\subseteq[d]\),

\[
 \operatorname{Pf}\!\left((v^2J_d)[A]\right)=v^{|A|}.
 \tag{5.2}
\]

**Proof.** The empty Pfaffian is one. Expanding the first row of \(J_{2r}\),
every remaining principal matrix is \(J_{2r-2}\), and every first-row entry is
one. Therefore

\[
 \operatorname{Pf}(J_{2r})=
 \left(\sum_{j=2}^{2r}(-1)^j\right)
 \operatorname{Pf}(J_{2r-2})
 =\operatorname{Pf}(J_{2r-2}),
\]

because the alternating sum has one more positive term than negative terms.
This proves (5.1) inductively. A principal submatrix \(J_d[A]\), in the
induced order, is \(J_{|A|}\); scaling a \(2r\times2r\) skew matrix by \(v^2\)
scales its Pfaffian by \((v^2)^r\). This gives (5.2). ∎

**[THEOREM — the degree-six local matchgate].** The signature (2.1) is the
pure-even sub-Pfaffian tensor of \(v^2J_6\), with empty entry one. Hence it is
a pure-even matchgate signature in the stated sub-Pfaffian convention.

The matrix stored in the artifact is explicitly

\[
 J_6=\begin{pmatrix}
0&1&1&1&1&1\\
-1&0&1&1&1&1\\
-1&-1&0&1&1&1\\
-1&-1&-1&0&1&1\\
-1&-1&-1&-1&0&1\\
-1&-1&-1&-1&-1&0
\end{pmatrix}.
\]

## 6. Exhaustive matchgate identities at \(d=6\)

For \(\alpha,\beta\in\{0,1\}^6\), let
\(p_1<\cdots<p_\ell\) be the positions where they differ. The matchgate
identities in the convention used by the producer are

\[
 \sum_{i=1}^{\ell}(-1)^i
 g(\alpha\mathbin\oplus e_{p_i})
 g(\beta\mathbin\oplus e_{p_i})=0.
 \tag{6.1}
\]

They can also be seen directly here. If \(\alpha,\beta\) are not both odd,
every product in (6.1) contains an odd entry of \(g\), so it vanishes. If both
are odd, then \(\ell\) is even. At a differing bit, the two Hamming-weight
changes cancel, so every nonzero product equals
\(v^{|\alpha|+|\beta|}\). The alternating sum of \(\ell\) equal terms is zero.

**[COMPUTATION].** At each of \(v=1/3\) and \(v=1/2\), all
\(64^2=4096\) ordered pairs are replayed with `Fraction` arithmetic. There are
992 algebraically active ordered identities: 480 with two terms, 480 with four
terms, and 32 with six terms. Every residual is exactly zero. The producer
uses recursive first-row Pfaffians; the verifier independently sums perfect
matchings with signs determined by chord crossings.

## 7. Ordinary real CP rank and every flattening

Let \(a=(1,v)\) and \(b=(1,-v)\). Entrywise,

\[
 g=\frac12\left(a^{\otimes6}+b^{\otimes6}\right),
 \tag{7.1}
\]

so the ordinary real CP rank is at most two.

Take any nonempty proper bipartition \(S\mid S^c\). Its flattening is

\[
 M_{S\mid S^c}=\frac12
 \left(a_Sa_{S^c}^{\mathsf T}+b_Sb_{S^c}^{\mathsf T}\right),
\]

where \(a_S(x)=v^{|x|}\) and \(b_S(x)=(-v)^{|x|}\). On every nonempty side,
these two vectors are independent: at the zero state they both equal one,
while at any one-bit state they equal \(v\) and \(-v\). Since \(v>0\), the
corresponding \(2\times2\) determinant is nonzero. Both side matrices have
column rank two, so the flattening has rank two. Thus:

**[THEOREM].** The ordinary real CP rank of \(g\) is exactly two, and every one
of the 62 ordered nontrivial bipartition flattenings has rank two.

For the first-three versus last-three split, identify the row and column
three-bit spaces. Let

\[
 u_{\rm e}(x)=\begin{cases}v^{|x|},&|x|\text{ even},\\0,&|x|\text{ odd},\end{cases}
\quad
 u_{\rm o}(x)=\begin{cases}0,&|x|\text{ even},\\v^{|x|},&|x|\text{ odd}.
\end{cases}
\]

Then

\[
 M_{3\mid3}=u_{\rm e}u_{\rm e}^{\mathsf T}
             +u_{\rm o}u_{\rm o}^{\mathsf T}.
\]

The vectors have disjoint parity support. Hence this operator is positive
semidefinite, has one rank-one even block and one rank-one odd block, and its
two nonzero eigenvalues are

\[
 \lambda_{\rm e}=1+3v^4,
 \qquad
 \lambda_{\rm o}=3v^2+v^6.
\]

The other six eigenvalues are zero.

**[COMPUTATION].** A self-contained `Fraction` Gauss-Jordan routine in the
producer checks all 62 flattenings at both rational audit points. The verifier
does not reuse it: it clears row denominators and performs integer
fraction-free cross-elimination.

## 8. Exact nonnegative all-leg CP rank

A nonzero nonnegative rank-one atom

\[
 T(x_1,\ldots,x_6)=\prod_{i=1}^6 f_i(x_i),\qquad f_i\ge0,
\]

has support \(R_1\times\cdots\times R_6\), with each nonempty
\(R_i\subseteq\{0,1\}\). There are three choices per leg and therefore
\(3^6=729\) nonempty Cartesian support rectangles.

**[LEMMA — rectangle barrier].** If such a rectangle is not a singleton, at
least one coordinate is free. Toggling that coordinate pairs each point with
a point of the opposite parity. Thus every nonsingleton rectangle contains
both even and odd states.

In a nonnegative sum representing \(g\), every atom must vanish at all 32 odd
states because the target is zero there and cancellation is impossible. The
lemma forces every nonzero atom to have support equal to one even singleton.
Since \(g\) is strictly positive at all 32 even states, at least 32 atoms are
necessary. Conversely,

\[
 g=\sum_{\substack{x\in\{0,1\}^6\\ |x|\text{ even}}}v^{|x|}
       \bigotimes_{i=1}^6\delta_{x_i}
 \tag{8.1}
\]

uses exactly 32 nonnegative atoms.

**[THEOREM].** The nonnegative all-leg CP rank of the Walsh-basis tensor \(g\)
is exactly 32.

**[COMPUTATION].** All 729 rectangles are enumerated. By their number of free
coordinates, the counts are

| free coordinates | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| rectangles | 64 | 192 | 240 | 160 | 60 | 12 | 1 |

Of the 64 singletons, 32 are even and 32 odd. Every one of the remaining 665
rectangles contains both parities. The verifier decodes the integers
\(0,\ldots,728\) in base three and rebuilds this inventory independently.

## 9. Why the positive spin-basis rank is only undoing decimation

Before the Walsh basis change,

\[
 \frac{W_6(\sigma)}{\cosh(K)^6}
 =\prod_{i=1}^6(1+v\sigma_i)
  +\prod_{i=1}^6(1-v\sigma_i).
 \tag{9.1}
\]

Every one-leg factor is strictly positive for \(0<v<1\). The two products are
the two values \(\tau=+1\) and \(\tau=-1\) of the original center spin. Their
factor vectors are not proportional when \(v>0\), so a nontrivial flattening
has rank two; the positive spin-basis CP rank is exactly two.

This does not contradict Section 8: nonnegative tensor rank is not invariant
under a signed Walsh basis change. More importantly, (9.1) introduces exactly
the hidden variable that decimation removed. It is evidence against, not for,
an arbitrary-hidden-auxiliary no-go.

## 10. Finite computation, provenance, and scope

The finite bounds stored in the artifact are:

- 4096 Walsh summands per rational audit point;
- 32 even principal Pfaffians per point, with at most 15 perfect-matching
  terms in one \(6\times6\) Pfaffian;
- 4096 ordered matchgate identities per point;
- 62 ordered nontrivial flattenings per point, with largest square core
  \(8\times8\);
- 729 Cartesian rectangles, with at most 64 points in one rectangle;
- at most 64 edge subsets, 16 local-star assignments, and 32 spin states in a
  cycle-space control;
- 84 total spin states and six nonempty even-sector tables in the
  \(d=2,4,6\) log-Walsh controls.

All mathematical values are integers or `Fraction`s. No fitted value, float
claim, benchmark, or critical coupling is used. The producer records current
SHA-256 hashes of itself, the verifier, and
`experiments/e234_star_decimation.py`; the verifier gates all three hashes and
does not import the producer. The producer writes JSON only after every check
passes. Both programs enforce a 2 GiB RSS wall. On Darwin they measure
`mach_task_basic_info.resident_size_max` through `task_info` flavor 20, avoiding
an inherited `ru_maxrss`; on Linux they use `ru_maxrss*1024`.

**[THEOREM — proved scope].** This front proves the local degree-six star
signature, its explicit pure-even Pfaffian certificate, all its ordinary
flattening ranks, its exact Walsh-basis nonnegative CP rank, the finite
bipartite cycle-space factorization, the positive spin-basis reconstruction,
and the all-even-degree log-Walsh sign law (4.4).

**[UNRESOLVED — not proved].** It does not prove:

- a no-go for arbitrary hidden auxiliary variables or unrestricted tensor
  networks;
- that a nonplanar or three-dimensional contraction of the local matchgate
  signatures is one global Pfaffian;
- global Pfaffian solvability or nonsolvability of the three-dimensional Ising
  partition function;
- integrability, a critical coupling, the thermodynamic free energy, or any
  critical exponent.

A local factor being Pfaffian does not erase the ordering, crossing, and
planarity obligations of a global matchgate contraction. Equation (3.1) is an
honest finite factor-graph identity, not a global solution theorem.

## 11. Sources for terminology and the one external product

The local algebra and all finite certificates above are proved here. The
following sources fix the standard terminology and identity conventions:

1. L. G. Valiant, “Quantum circuits that can be simulated classically in
   polynomial time,” *SIAM Journal on Computing* **31** (2002), 1229–1254,
   DOI: [10.1137/S0097539700377025](https://doi.org/10.1137/S0097539700377025).
   This is primary matchgate literature.
2. J.-Y. Cai and A. Gorenstein, “Matchgates Revisited,” arXiv:1303.6729,
   especially Theorem 1 for (6.1) and the necessity/sufficiency discussion:
   [https://arxiv.org/abs/1303.6729](https://arxiv.org/abs/1303.6729).
3. G. D. Forney, Jr. and P. O. Vontobel, “Partition Functions of Normal Factor
   Graphs,” arXiv:1102.0316 / 2011 Information Theory and Applications
   Workshop: [https://arxiv.org/abs/1102.0316](https://arxiv.org/abs/1102.0316).
   This is the NFG sum-of-products convention used in Section 3.
4. NIST Digital Library of Mathematical Functions, Eq. 4.36.E2, the product
   (4.2): [https://dlmf.nist.gov/4.36.E2](https://dlmf.nist.gov/4.36.E2).
