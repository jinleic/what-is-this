# Exact sixth-order interlayer reduction to 2D connected correlations

Artifact: `results/interlayer/c6_closed.json`  
Producer: `experiments/e99_c6_closed.py`  
Clean-room verifier: `tests/test_c6_closed.py`

All coefficients in this note and both certificates are `int` or
`fractions.Fraction`; no floating-point calculation enters the result.

## 1. Statement

Let \(\langle\cdot\rangle_2\) denote the zero-field square-lattice Ising
layer state, interpreted coefficientwise as a high-temperature formal series.
For six labeled in-plane sites \(x_1,\ldots,x_6\), set \(x_1=0\), and write

\[
 M_I=\left\langle\prod_{i\in I}\sigma_{x_i}\right\rangle_2.
\]

Repeated physical positions are permitted: their spins cancel in the product,
or equivalently the associated site masks are XORed.  Define the connected 2D
correlations

\[
\begin{aligned}
 G_{ij}&=M_{ij},\\
 U_{ijkl}&=M_{ijkl}-G_{ij}G_{kl}-G_{ik}G_{jl}-G_{il}G_{jk},\\
 W_{123456}&=M_{123456}
 -\sum_{\{ij\}\sqcup R=\{1,\ldots,6\}}G_{ij}U_R
 -\sum_{\{ij\}\sqcup\{kl\}\sqcup\{mn\}}G_{ij}G_{kl}G_{mn}.
\end{aligned}
\]

The first sum has 15 labeled \(2+4\) splits and the second has 15 pairings.
For every displayed term below, let \(\sum'\) mean
\(\sum_{x_2,\ldots,x_6\in\mathbb Z^2}\), coefficientwise in \(v=\tanh K\).

[THEOREM] The direct-coupling sixth coefficient of the stacked independent
2D layers is the following finite linear combination of 2D connected
correlation lattice sums:

\[
\begin{aligned}
 c_{6,q}=\sum'\Big[&
 \frac1{720}W_{123456}^2
 +\frac1{72}G_{12}U_{3456}W_{123456}
 +\frac5{72}G_{23}U_{1456}W_{123456}\\
 &+\frac1{36}G_{12}G_{13}U_{2456}U_{3456}
 +\frac5{36}G_{12}G_{23}U_{1456}U_{3456}
 +\frac16G_{12}G_{34}G_{56}W_{123456}\\
 &+\frac13G_{12}G_{34}U_{1256}U_{3456}
 +\frac13G_{23}G_{24}U_{1356}U_{1456}
 +\frac16G_{23}G_{45}U_{1236}U_{1456}\\
 &+\frac{11}{12}G_{12}G_{13}G_{24}G_{56}U_{3456}
 +\frac{11}{12}G_{12}G_{23}G_{34}G_{56}U_{1456}\\
 &+\frac76G_{12}G_{34}G_{35}G_{46}U_{1256}
 +\frac53G_{12}G_{13}G_{24}G_{35}G_{46}G_{56}
 \Big].
\tag{1}
\end{aligned}
\]

The notation in (1) is canonical only under relabeling of the five free
slots.  The coefficient of each representative already includes its full
\(S_5\) orbit, because the primed sum is invariant under that relabeling.
Thus (1) is an identity of formal series, not a claim of absolute convergence
of the displayed infinite sums at a specified temperature.

## 2. Partition-block derivation

Put \(B_g(x)=\sigma_{g,x}\sigma_{g+1,x}\), the vertical bond at in-plane
site \(x\) across gap \(g\).  Connected vertical support and independent
layer spin-flip symmetry leave just three sixth-order gap profiles.  With
\(Y_i=B_0(x_i)\), \(Z_i=B_1(x_i)\), and \(T_i=B_2(x_i)\), their Taylor
placement factors give

\[
\begin{aligned}
 c_{6,q}=\sum'\bigg[&\frac{1}{6!}\kappa(Y_1,Y_2,Y_3,Y_4,Y_5,Y_6)\\
 &+\frac{2}{4!2!}\kappa(Y_1,Y_2,Y_3,Y_4,Z_5,Z_6)\\
 &+\frac{1}{(2!)^3}\kappa(Y_1,Y_2,Z_3,Z_4,T_5,T_6)\bigg].
\tag{2}
\end{aligned}
\]

The factor two in the middle line combines reflected \(4+2\) and \(2+4\)
profiles.  A skipped gap splits the layer sets into independent families, and
an odd incident count in any layer is killed by that layer's spin flip.

For a labeled set partition \(\pi\), the cumulant expansion is

\[
 \kappa(A_1,\ldots,A_6)=
 \sum_{\pi}\mu(\pi)\prod_{C\in\pi}
 \left\langle\prod_{i\in C}A_i\right\rangle,
 \qquad
 \mu(\pi)=(-1)^{|\pi|-1}(|\pi|-1)! .
\tag{3}
\]

There are exactly \(B_6=203\) partitions.  For a block \(C\), factorization
over independent 2D layers makes its moment respectively

\[
 M_C^2,\qquad
 M_{C_0}M_{C_1}M_{C_0\cup C_1},\qquad
 M_{C_0}M_{C_2}M_{C_0\cup C_1}M_{C_1\cup C_2}
\tag{4}
\]

for the one-, two-, and three-gap profiles.  Here \(C_j\) is the subset of
slots on gap \(j\); when positions coincide, each \(M\) in (4) is evaluated
using the XOR site mask as specified above.  Formula (4) is therefore valid
without excluding coincident indices.

After zero layer moments are removed, the numbers of admissible partitions in
the three profiles are \(31,11,5\).  Before connected-correlation collection,
the exact expansion has 376 monomials.  Substituting the displayed
moment-to-\(G/U/W\) identities and collecting under \(S_5\) gives precisely
the 13 terms in (1).  The profile Möbius-weight checks are

\[
\begin{array}{c|rrr}
\text{gap multiplicities}&6&4+2&2+2+2\\ \hline
6&1&-15&30\\
4+2&1&-7&6\\
2+4&1&-7&6\\
2+2+2&1&-3&2
\end{array}
\]

and the resulting normal-form SHA-256 is
`44510d76fcea9c6cc60c5e6f9858b99de150022fe7df9162f06e033a76374eb3`.
This is an exact finite partition-algebra certificate of (1).

## 3. Exact formal-series certificate through \(v^{12}\)

[COMPUTATION] For every open \(a\times b\) rectangle, the producer computes
all layer moments as \(m_S=P_S/P_\varnothing\), where \(P_S\) is the exact
integer even-subgraph boundary polynomial.  It forms the raw vertical subset
columns from (4), applies (3), performs vertical Möbius inversion through
heights 2, 3, and 4, and then planar finite-lattice inversion over all 28
rectangles satisfying \((a-1)+(b-1)\leq6\).  This independently re-expands
the raw partition form underlying (1), using only exact formal arithmetic.

The three connected vertical-profile rows are

\[
\begin{array}{c|rrrrrrr}
 &v^0&v^2&v^4&v^6&v^8&v^{10}&v^{12}\\ \hline
6&1/45&34/45&202/15&7406/45&74876/45&684746/45&2043622/15\\
4+2,2+4&0&-8/3&-332/3&-6184/3&-25312&-208248&-1845892/3\\
2+2+2&0&2&158&5102&114696&2038010&30877538
\end{array}
\]

They sum coefficientwise to

\[
\boxed{
 c_{6,q}=\frac1{45}+\frac4{45}v^2+\frac{304}{5}v^4
 +\frac{144236}{45}v^6+\frac{4097156}{45}v^8
 +\frac{83024036}{45}v^{10}+\frac{455977232}{15}v^{12}
 +O(v^{14}).}
\tag{5}
\]

The exact change of interlayer variable is

\[
 c_{6,w}^{\rm total}=c_{6,q}+\frac43c_{4,q}+\frac{23}{45}c_{2,q}.
\tag{6}
\]

After subtracting the separate vertical \([w^6]\log\cosh K_z=1/6\)
prefactor, (6) gives

\[
 c_{6,w}^{\rm residual}=2v^2+138v^4+\frac{13682}{3}v^6+109718v^8
 +2061698v^{10}+32654478v^{12}+O(v^{14}),
\]

exactly matching every stored coefficient in `results/interlayer/c6_series.json`.

## 4. Standalone verification

`tests/test_c6_closed.py` imports neither producer experiment nor producer
helper.  It separately:

1. enumerates all 203 partitions and reconstructs the 13-term normal form;
2. exhaustively evaluates both the raw profile expression (2) and (1) on an
   open \(2\times2\) layer over every \(4^5\) anchored position tuple,
   coefficientwise through \(v^{12}\);
3. rebuilds the boundary-mask transfer, vertical cumulants, and planar FLM,
   then matches every direct-\(q\), profile, and converted-\(w\) coefficient
   through \(v^{12}\).

Reproduce with:

```sh
timeout 1800 .venv/bin/python experiments/e99_c6_closed.py
timeout 1800 .venv/bin/python tests/test_c6_closed.py
```

Both commands print `PASS` only after the exact checks described above.

## 5. Scope

[UNRESOLVED] Equation (1) reduces the sixth interlayer cumulant to standard
2D two-, four-, and six-point connected correlation objects.  The latter have
Onsager/Kaufman fermionic/Pfaffian evaluation in principle, but this work does
not provide a scalar closed evaluation of the resulting infinite 2D lattice
sums.  It makes no analytic-continuation, convergence-radius, critical-point,
or 3D-solution claim.
