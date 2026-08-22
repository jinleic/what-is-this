# Sixth-order interlayer coefficient from exact layer cumulants

## Scope and conventions

[LEMMA] Let
\[
v=\tanh K,\qquad q=K_z,\qquad w=\tanh K_z,
\]
and let
\[
X_{\ell i}=\sigma_{\ell i}\sigma_{\ell+1,i},\qquad
V=\sum_{\ell,i}X_{\ell i}.
\]
The cumulant expansion is an expansion in the direct coupling `q`:
\[
\frac1{|\Lambda|}\log\langle e^{qV}\rangle_0
 =\sum_{n\geq 1}\frac{q^n}{n!|\Lambda|}\kappa_n(V),
\]
where the expectation is over independent zero-field square-lattice layers. The high-temperature graph variable is instead `w`. The two variables are not identified in this calculation.

[LEMMA] Through sixth order,
\[
q=\operatorname{atanh}w=w+\frac{w^3}{3}+\frac{w^5}{5}+O(w^7),
\]
so
\[
q^2=w^2+\frac23w^4+\frac{23}{45}w^6+O(w^8),\qquad
q^4=w^4+\frac43w^6+O(w^8),\qquad q^6=w^6+O(w^8).
\]
Consequently,
\[
\boxed{c_{6,w}^{\rm tot}=c_{6,q}+\frac43c_{4,q}+\frac{23}{45}c_{2,q}.}
\]
The vertical-bond high-temperature prefactor is
\[
\log\cosh(\operatorname{atanh}w)
=-\frac12\log(1-w^2)
=\frac12w^2+\frac14w^4+\frac16w^6+O(w^8),
\]
and hence `c6_w_res = c6_w_tot - 1/6` at `v^0`.

[LEMMA] An open stack is invariant under flipping every other layer, which sends every `X_(ell i)` to its negative while leaving all in-plane energies unchanged. Thus the interlayer free energy is even in `q`, equivalently even in `w`. In particular `c3=c5=0` exactly, not merely through the computed in-plane order.

## Complete sixth-cumulant enumeration

[LEMMA] For six labeled random variables the complete joint-cumulant formula is
\[
\boxed{
\kappa(Y_1,\ldots,Y_6)=
\sum_{\pi\in\Pi_6}(-1)^{|\pi|-1}(|\pi|-1)!
\prod_{B\in\pi}\left\langle\prod_{a\in B}Y_a\right\rangle .}
\]
There are `B_6=203` set partitions. The experiment generates all of them recursively and only then groups them by block-size profile. No sixth-order subtraction coefficient is hand-entered.

[LEMMA] For the symmetric vertical-bond sum, odd raw moments vanish. Grouping the 203 generated partitions leaves
\[
\boxed{\kappa_6(V)=\mu_6-15\mu_4\mu_2+30\mu_2^3,}
\qquad \mu_j=\langle V^j\rangle.
\]
The three surviving profiles contain respectively `1`, `15`, and `15` set partitions; their total Möbius weights are `1`, `-15`, and `30`. The artifact also records all eight profiles killed by an odd block, so the reduction can be audited against all 203 inputs.

[LEMMA] The layer-moment form used before any summation over sites is as follows. Assign vertical bond `a` to gap `g_a` and site `i_a`. For a block `B` of a set partition, define the spin mask in layer `ell` by symmetric difference,
\[
S_\ell(B)=\mathop{\triangle}_{\substack{a\in B\\g_a=\ell-1\text{ or }g_a=\ell}}\{i_a\},
\]
and write `M(S)=<prod_(i in S) sigma_i>_2D`. Independence of layers gives
\[
\left\langle\prod_{a\in B}X_{g_a i_a}\right\rangle_0
=\prod_\ell M(S_\ell(B)).
\]
Therefore the exact fixed-gap, fixed-site sixth cumulant is
\[
\boxed{
\sum_{\pi\in\Pi_6}(-1)^{|\pi|-1}(|\pi|-1)!
\prod_{B\in\pi}\prod_\ell M(S_\ell(B)).}
\]
This expression also handles repeated sites correctly because repeated insertions cancel in the symmetric difference.

### Exhaustive connected layer-gap profiles

[LEMMA] For any nonzero block moment, every layer contains an even number of spin insertions. Starting at an outer layer and moving inward shows that each gap used by that block is crossed an even number of times. A connected sixth-order term therefore uses at most three layer gaps. Cumulant connectedness eliminates a skipped gap. After anchoring the lowest used gap, the exhaustive profiles are
\[
(6),\qquad(4,2),\qquad(2,4),\qquad(2,2,2).
\]
They have exact vertical extents two, three, three, and four layers. There are `binom(6,2)=15` labeled gap assignments for each of `(4,2)` and `(2,4)`, and `6!/(2!2!2!)=90` assignments for `(2,2,2)`.

[LEMMA] The same 203-partition generator was specialized to each fixed gap assignment. A partition survives outer-layer parity only when every one of its blocks crosses every gap an even number of times. The resulting audit table is

| ordered gap counts | admissible `6` blocks | admissible `4+2` partitions | admissible `2+2+2` partitions | summed Möbius weights |
|---|---:|---:|---:|---|
| `6` | 1 | 15 | 15 | `1, -15, 30` |
| `4+2` | 1 | 7 | 3 | `1, -7, 6` |
| `2+4` | 1 | 7 | 3 | `1, -7, 6` |
| `2+2+2` | 1 | 3 | 1 | `1, -3, 2` |

The counts do not identify different site-moment factors; the computation retains the labeled blocks and evaluates their layer masks before summation. They are an omission certificate for every adjacent-gap subtraction term.

[LEMMA] A profile supported on two separated gap intervals factors into two independent layer families. The generated 203-term polynomial was evaluated with a distinct symbolic moment for every subset in each family and expands identically to zero. Thus disconnected gap profiles are removed by the cumulant rather than discarded numerically.

## Exact two-dimensional moment engine

[LEMMA] For an open square rectangle `R`, let
\[
P_S(v)=\sum_{E\subseteq E(R):\,\partial E=S}v^{|E|}.
\]
The exact layer moment is
\[
M(S)=\left\langle\prod_{i\in S}\sigma_i\right\rangle_R
=\frac{P_S(v)}{P_\varnothing(v)}.
\]
The experiment enumerates every `P_S` with integers through `v^12`; every subsequent division and cumulant operation uses `Fraction` arithmetic.

[LEMMA] For a height-`h` open stack, introduce a boundary mask `S_j` on each of its `h-1` gaps. Summing the independent-layer moments gives the normalized vertical-subset polynomial
\[
\mathcal Z_h(w)=\frac{1}{P_\varnothing(v)^h}
\sum_{S_1,\ldots,S_{h-1}}
 w^{\sum_j|S_j|}
 P_{S_1}P_{S_1\triangle S_2}\cdots
 P_{S_{h-2}\triangle S_{h-1}}P_{S_{h-1}}.
\]
This is evaluated through vertical degree six by a transfer over the last mask. It includes every same-gap and adjacent-gap mask sequence.

[LEMMA] If `A=|R|` and `B=A(h-1)`, then
\[
\left\langle e^{qV_h}\right\rangle_0
=\cosh(q)^B\mathcal Z_h(\tanh q).
\]
Thus each raw moment is obtained exactly as
\[
\mu_n(h)=n![q^n]\,\cosh(q)^B\mathcal Z_h(\tanh q),
\qquad 1\le n\le6.
\]
The generated set-partition formula then supplies `kappa_n(h)/n!`.

[LEMMA] Vertical interval Möbius inversion isolates exact height:
\[
W_n(h)=\frac{\kappa_n(V_h)}{n!}
-\sum_{r<h}(h-r+1)W_n(r).
\]
It yields height two only for `n=2`, heights two and three for `n=4`, and heights two, three, and four for `n=6`. The computed forbidden weights `W_2(3)`, `W_2(4)`, and `W_4(4)` vanish coefficientwise through `v^12`, agreeing with the cut-parity lemma.

[LEMMA] Planar rectangular Möbius inversion removes all terms embeddable in a proper subrectangle. A connected even graph spanning an `a` by `b` in-plane box crosses every separating coordinate cut a positive even number of times. It therefore has at least
\[
2((a-1)+(b-1))
\]
in-plane edges. Through `v^12`, all required rectangles obey `(a-1)+(b-1)<=6`. The computation uses all 28 oriented rectangles in that set. Adding every rectangle with one extra span unit changes no component through `v^8`.

## Independent anisotropic spin-DOS route

[COMPUTATION] The independent route does not import the layer-moment experiment. For each open three-dimensional box of height one through four and in-plane span at most three, the repository's exact joint spin density-of-states engine is transformed directly to the anisotropic even-subgraph polynomial
\[
P(v,w)=P_0(v)+P_2(v)w^2+P_4(v)w^4+P_6(v)w^6+O(w^8).
\]
It then extracts
\[
\begin{aligned}
[w^2]\log P&=r_2,\\
[w^4]\log P&=r_4-\frac12r_2^2,\\
[w^6]\log P&=r_6-r_2r_4+\frac13r_2^3,
\end{aligned}
\qquad r_j=P_j/P_0,
\]
and performs an independent three-dimensional box Möbius inversion.

[COMPUTATION] This spin-DOS route agrees coefficientwise through `v^6` with the layer-moment/set-partition route for `c2`, `c4`, and `c6`. Its exact height-two, height-three, and height-four `w^6` weights separately agree with the corresponding converted cumulant components. Every enumerated `w^1`, `w^3`, and `w^5` column is zero.

## Exact coefficients

[COMPUTATION] The generic cumulant engine reproduces the established lower orders:
\[
\begin{aligned}
c_{2,q}(v)&=\frac12+2v^2+18v^4+118v^6+778v^8+4978v^{10}+31398v^{12}+O(v^{14}),\\
c_{4,q}(v)&=-\frac1{12}+\frac23v^2+51v^4+\frac{2914}{3}v^6
+\frac{41113}{3}v^8+\frac{481894}{3}v^{10}+1679961v^{12}+O(v^{14}).
\end{aligned}
\]

[COMPUTATION] The same-gap, exact-height-two contribution to the direct-coupling sixth coefficient is
\[
\begin{aligned}
c_{6,q}^{(2)}(v)={}&\frac1{45}+\frac{34}{45}v^2+\frac{202}{15}v^4
+\frac{7406}{45}v^6+\frac{74876}{45}v^8\\
&+\frac{684746}{45}v^{10}+\frac{2043622}{15}v^{12}+O(v^{14}).
\end{aligned}
\]

[COMPUTATION] The two-adjacent-gap, exact-height-three contribution, including both ordered profiles `(4,2)` and `(2,4)`, is
\[
\begin{aligned}
c_{6,q}^{(3)}(v)={}&-\frac83v^2-\frac{332}{3}v^4-\frac{6184}{3}v^6
-25312v^8\\
&-208248v^{10}-\frac{1845892}{3}v^{12}+O(v^{14}).
\end{aligned}
\]

[COMPUTATION] The three-adjacent-gap, exact-height-four contribution is
\[
\begin{aligned}
c_{6,q}^{(4)}(v)={}&2v^2+158v^4+5102v^6+114696v^8\\
&+2038010v^{10}+30877538v^{12}+O(v^{14}).
\end{aligned}
\]

[COMPUTATION] Summing all exact vertical extents gives the new direct-coupling prefix
\[
\boxed{\begin{aligned}
c_{6,q}(v)={}&\frac1{45}+\frac4{45}v^2+\frac{304}{5}v^4
+\frac{144236}{45}v^6+\frac{4097156}{45}v^8\\
&+\frac{83024036}{45}v^{10}+\frac{455977232}{15}v^{12}+O(v^{14}).
\end{aligned}}
\]

[COMPUTATION] After the exact `q`-to-`w` conversion, the total coefficient including the single-bond prefactor is
\[
\boxed{\begin{aligned}
c_{6,w}^{\rm tot}(v)={}&\frac16+2v^2+138v^4+\frac{13682}{3}v^6
+109718v^8\\
&+2061698v^{10}+32654478v^{12}+O(v^{14}).
\end{aligned}}
\]
Therefore the residual after separating `log cosh(K_z)` is
\[
\boxed{c_{6,w}^{\rm res}(v)=
2v^2+138v^4+\frac{13682}{3}v^6+109718v^8
+2061698v^{10}+32654478v^{12}+O(v^{14}).}
\]

[COMPUTATION] The independent spin-DOS route obtains
\[
c_{6,w}^{\rm res}(v)=2v^2+138v^4+\frac{13682}{3}v^6+O(v^8)
\]
exactly. The standalone test recomputes this route without importing the experiment.

## Verification and limitation

[COMPUTATION] `experiments/e64_interlayer_c6.py` prints `PASS` only after all embedded exact checks pass and writes `results/interlayer/c6_series.json` with `provenance/data/checks`. `tests/test_interlayer_c6.py` independently recomputes the direct spin-DOS finite-lattice route through `v^6` and prints `PASS test_interlayer_c6`.

[UNRESOLVED] These are exact coefficients of a finite high-temperature expansion around `K=K_z=0`. They do not establish an all-orders formula, determine a convergence radius, justify substituting the benchmark critical coupling, or solve the three-dimensional Ising model.
