# Two exact finite obstructions beyond the infrared endpoint

Artifacts: `experiments/e175_upper_endpoint4_paired.py`,
`experiments/e176_upper_endpoint4_current.py`,
`experiments/e177_upper_endpoint4.py`,
`results/bounds/upper_endpoint4.json`, and
`tests/test_upper_endpoint4.py`.

## 0. Outcome and exact scope

**[UNRESOLVED]** No improved upper endpoint is obtained.  The certified upper
endpoint remains

\[
 K_c\le \bar K:=I_3/2
 =0.2527310098586630030260020266135701299926\ldots .
\]

**[COMPUTATION]** Both finite certificate classes are challenged at the
benchmark-independent rational high-temperature point

\[
 v_\dagger=\frac6{25},\qquad
 K_\dagger=\operatorname {atanh}(v_\dagger)
 \in[0.244774112659352892960068,
      0.244774112659352892960069].                    \tag{0.1}
\]

**[COMPUTATION]** Exact positive-series remainders and the inherited rational
lower endpoint for `I3` prove

\[
 K_\dagger<\bar K,\qquad
 \bar K-K_\dagger>
 0.007956897199310110065933.                           \tag{0.2}
\]

The exact rational endpoints of (0.1)--(0.2), rather than the displayed
decimals, are stored in the artifact.

**[THEOREM]** In the finite class `PM4-aggregate-one-moment`, the exact SDP
optimum is

\[
 \sup s=1,\qquad \eta_*:=1-\sup s=0.                  \tag{0.3}
\]

Thus this smallest aggregate genuinely four-point relaxation supplies no
uniform magnetisation floor.  Its pre-fixed rational challenge
`eta=1/1024` is infeasible by the exact gap `1/1024`.

**[THEOREM]** In the finite class
`RC2-single-edge-conditional-domination`, the exact LP optimum at
`v=6/25` is

\[
 p_*=v_\dagger^2=\frac{36}{625}.                       \tag{0.4}
\]

The self-contained planar-percolation certificate below requires
`p>=5/6`; the exact infeasibility gap is

\[
 \frac56-p_*=\frac{2909}{3750}.                        \tag{0.5}
\]

**[UNRESOLVED]** Equations (0.3) and (0.5) are obstructions only to the two
named finite classes.  They do not exclude mode-resolved four-point or DLR
constraints, multi-edge/block random-current certificates, source-dependent
current inequalities, or any other extension outside those definitions.

## 1. Inherited walls and the challenge point

**[EXTERNAL]** `proofs/upper_infrared.md` proves that the two-point
reflection-positivity/GKS/infrared class is exactly optimal at `I3/2`;
`proofs/mag_floor.md` proves that its finite-torus LP floor is at most
`5168/(525L)` for the stated range; and `proofs/upper_beyond.md` proves that
the plain Peierls head-plus-tail class cannot reach the incumbent.  No part of
those artifacts is modified here.

**[LEMMA]** For `0<v<1` and an integer `m>=1`,

\[
 \sum_{n=0}^{m-1}\frac{v^{2n+1}}{2n+1}
 \le \operatorname {atanh}v
 \le \sum_{n=0}^{m-1}\frac{v^{2n+1}}{2n+1}
 +\frac{v^{2m+1}}{(2m+1)(1-v^2)} .                    \tag{1.1}
\]

**[LEMMA]** The lower inequality in (1.1) drops a positive tail, while in the
tail `1/(2n+1)<=1/(2m+1)` and the remaining powers form a geometric series.
Thus (1.1) uses exact rational arithmetic at rational `v`.

**[COMPUTATION]** Taking `m=20`, `v=6/25`, and comparing `2*K_hi` with the
stored rational `I3_lo` proves (0.1)--(0.2).  The independent verifier uses
`m=13` and also checks the coarser different bound

\[
 \operatorname {atanh}v\le v+\frac{v^3}{3(1-v^2)},
 \qquad 2\left(v_\dagger+
 \frac{v_\dagger^3}{3(1-v_\dagger^2)}\right)
 =\frac{7212}{14725}<I_{3,\mathrm{lo}}.                \tag{1.2}
\]

**[COMPUTATION]** At this target the inherited paired-defect threshold is

\[
 \delta_\dagger=\frac{I_3}{2K_\dagger}-1
 =0.03250710262152419\ldots>0,                         \tag{1.3}
\]

with an exact rational enclosure stored in the artifact.  Therefore the
thermodynamic-limit summed infrared cap on the nonzero-mode mass lies strictly above one and is
redundant against the exact bounded-spin cap used next.

## 2. Paired momentum: the exact smallest aggregate SDP

### 2.1 A genuinely four-point variable

**[LEMMA]** On an even torus of `N` sites, for every spin configuration define

\[
 R(\sigma)=\frac1{N^2}\sum_{k\ne0}|\widehat\sigma(k)|^2.
\]

Configuration-wise Parseval and `sigma_x^2=1` give

\[
 R(\sigma)=1-\left(\frac{M(\sigma)}N\right)^2,
 \qquad 0\le R\le1.                                   \tag{2.1}
\]

**[LEMMA]** Put

\[
 s=\mathbb E R=1-M_L^2,
 \qquad q=\mathbb E R^2
 =\frac1{N^4}\sum_{k,l\ne0}
   \mathbb E\!igl[|\widehat\sigma(k)|^2
                    |\widehat\sigma(l)|^2\bigr].      \tag{2.2}
\]

The variable `q` is a paired-momentum four-spin moment, not a two-point
repackaging.  Equations (2.1)--(2.2) imply exactly

\[
 \begin{pmatrix}1&s\\s&q\end{pmatrix}\succeq0,
 \qquad q\ge0,
 \qquad s-q\ge0,
 \qquad1-s\ge0.                                       \tag{2.3}
\]

The PSD row is `Var(R)>=0`, and `s-q=E[R-R^2]>=0`.

**[COMPUTATION]** This is the smallest scalar aggregate moment system that
contains a four-point variable: `s` is the quadratic moment, `q` its second
moment, and the `2x2` moment matrix is the first nontrivial PSD matrix.  The
class name `PM4-aggregate-one-moment` refers exactly to (2.3) together with
the redundant thermodynamic-limit summed infrared cap at `K_dagger`.

### 2.2 The finite certificate feasibility problem

**[THEOREM]** The primal SDP is the exact rational problem

\[
 \begin{array}{ll}
 \text{maximize}&s\\
 \text{subject to}&(2.3).
 \end{array}                                           \tag{2.4}
\]

A dual order-floor certificate with `eta>0` is a rational solution
`a,b,c>=0` and a rational PSD matrix
`Z=[[z00,z01],[z01,z11]]` of the coefficient identity

\[
 1-\eta-s
 =a q+b(s-q)+c(1-s)
  +\left\langle Z,
    \begin{pmatrix}1&s\\s&q\end{pmatrix}\right\rangle. \tag{2.5}
\]

Equivalently, the finite affine/semidefinite feasibility rows are

\[
 \begin{aligned}
 c+z_{00}&=1-\eta,\\
 b-c+2z_{01}&=-1,\\
 a-b+z_{11}&=0,\\
 a,b,c&\ge0,\qquad Z\succeq0.                          \tag{2.6}
 \end{aligned}
\]

**[THEOREM]** Any exact rational feasible point of (2.6) with `eta>0` proves
`s<=1-eta` for every physical torus point, because every term on the
right-hand side of (2.5) is nonnegative.  Hence
`M_L^2=1-s>=eta` uniformly in `L`; the thermodynamic implication already
proved in `proofs/upper_infrared.md` then gives

\[
 K_c\le K_\dagger<I_3/2.                              \tag{2.7}
\]

Thus rational feasibility is a rigorous endpoint certificate, not a numerical
suggestion.

### 2.3 Exact obstruction and measured gap

**[THEOREM]** The rational primal point

\[
 (s,q)=(1,1),\qquad
 \begin{pmatrix}1&1\\1&1\end{pmatrix}\succeq0          \tag{2.8}
\]

satisfies every row of (2.3).  The dual identity `1-s=1*(1-s)` gives the
matching upper bound `s<=1`.  Therefore the exact optimum of (2.4) is `1`,
and the exact optimum of (2.6) is `eta_*=0`.

**[THEOREM]** This paired obstruction also covers the whole incumbent range.
For every `K<=I3/2`, the thermodynamic-limit summed infrared row is
`s<=I3/(2K)` with right-hand side at least one, so it is redundant against
`s<=1`.  The same primal and dual points therefore give `eta_*(K)=0` for every
`K<=I3/2`.  The class cannot improve the incumbent at any coupling where an
improvement is needed.

**[COMPUTATION]** The witness (2.8) has a literal bounded-spin lift: on the
`2x2x2` checkerboard `sigma(x)=(-1)^(x_1+x_2+x_3)`, the zero Fourier power is
zero, the seven nonzero powers are `(0,0,0,0,0,0,64)`, and hence `R=R^2=1`.
An orbit mixture supplies the same aggregate moments if translation symmetry
is desired.

**[COMPUTATION]** For the fixed challenge `eta=1/1024`, (2.8) has
`s=1>1023/1024` by exactly `1/1024`.  Substitution in (2.5) makes its left
side `-1/1024` while the right side is nonnegative, an exact infeasibility
certificate.

**[UNRESOLVED]** This theorem does not cover separate mode variables,
mode-to-mode fourth cumulants, local DLR identities, or a four-point inequality
that excludes (2.8).  It closes only `PM4-aggregate-one-moment`.

## 3. Random currents: an exact one-edge domination LP

### 3.1 Switching weights with no transcendental arithmetic

**[LEMMA]** On a finite graph, let `n_1,n_2` be two independent sourceless
integer currents at coupling `K`, and put `v=tanh K`.  After dividing each
current's one-edge weight sum by `cosh K`, summing the counts while retaining
each colour parity and whether the union edge is present gives exactly five
states:

\[
\begin{array}{c|c|c}
\text{colour parities}&\text{union present?}&\text{weight}\\\hline
00&0&1-v^2\\
00&1&v^2\\
11&1&v^2\\
10&1&v\\
01&1&v.
\end{array}                                           \tag{3.1}
\]

**[LEMMA]** The first two weights are respectively
`sech(K)^2=1-v^2` and `1-sech(K)^2=v^2`; the third is
`tanh(K)^2=v^2`; and each mixed-parity weight is `tanh(K)=v`.  These identities
sum the full integer-current series exactly and leave rational weights at
rational `v`.

**[LEMMA]** Condition on all other coloured edge states.  The source constraints
fix the residual two-colour parity on the remaining edge.  In residual state
`00`, (3.1) gives conditional union-open probability `v^2`; in residual states
`01`, `10`, or `11`, the edge is certainly open.  Therefore every one-edge
conditional union-open probability is at least `v^2`.

**[LEMMA]** Revealing edges sequentially and coupling each reveal with a common
uniform random variable proves that the union support stochastically dominates
iid bond percolation with parameter `p=v^2`.  The same conclusion follows for
partial reveals because their conditional probability is an average of the
full-other-edge conditional probabilities.

**[LEMMA]** The switching involution along a deterministic union path from
`x` to `y` gives, on every finite graph,

\[
 \mathbb P^{\varnothing,\varnothing}(x\leftrightarrow y)
 =\langle\sigma_x\sigma_y\rangle^2.                   \tag{3.2}
\]

The path switch bijects two sourceless currents whose union connects `x,y`
with two currents carrying source set `{x,y}`; current weights are preserved.

### 3.2 A self-contained sufficient percolation threshold

**[LEMMA]** Iid bond percolation with `p>=5/6` on a square-lattice plane has
positive percolation probability.  At closed-edge probability `q=1/6`, the
number of length-`n` closed dual circuits surrounding a fixed vertex is at most
`4*n*3^(n-1)`.
Mark the first crossing of one of four coordinate half-axes: there are at most
`4n` marked starts, and after the first directed edge there are at most three
nonbacktracking choices per step.  Hence

\[
 \mathbb P(\text{a surrounding closed dual circuit})
 \le\frac43\sum_{n\ge4}n(3q)^n
 =\frac43\left[\frac{r^4(4-3r)}{(1-r)^2}\right]_{r=1/2}
 =\frac56<1.                                           \tag{3.3}
\]

Thus the probability that the vertex belongs to an infinite open cluster is at
least `1/6`.

**[LEMMA]** At `q=1/6`, the probability of a length-`n` closed dual self-avoiding
path from a fixed dual vertex is at most
`4*3^(n-1)q^n=(4/3)2^{-n}`, which tends to zero.  Hence there is almost surely
no infinite closed dual path.  Planarity then excludes two disjoint infinite
open clusters.  Product-measure FKG consequently gives, for any two vertices
`x,y` in the plane,

\[
 \mathbb P_{5/6}(x\leftrightarrow y)\ge(1/6)^2=1/36.   \tag{3.4}
\]

**[THEOREM]** If the two-current union at `K_dagger` admits a rational iid
one-edge domination parameter `p>=5/6`, then (3.2)--(3.4) give
`<sigma_x sigma_y> >= 1/6` uniformly along the plane.  Therefore the model is
ordered at `K_dagger` and (2.7) follows.  This implication uses the displayed
finite-current switching identity and the proved contour sum (3.3), not an
unproved asymptotic differential inequality.

### 3.3 The finite LP and its exact obstruction

**[THEOREM]** The class `RC2-single-edge-conditional-domination` is the exact
one-variable linear program

\[
 \begin{array}{ll}
 \text{maximize}&p\\
 \text{subject to}&p\ge0,\quad p\le v_\dagger^2,\quad
                    p\le1\quad\text{on each of rows }01,10,11.
 \end{array}                                           \tag{3.5}
\]

Adding the rational certificate row `p>=5/6` makes (3.5) a finite feasibility
problem whose rational feasible points rigorously imply (2.7).

**[THEOREM]** The primal point `p=36/625` is feasible, and the residual-`00`
row is a one-row exact dual upper certificate.  Thus

\[
 p_*=\frac{36}{625},\qquad
 \frac56-p_*=\frac{2909}{3750}>0,                     \tag{3.6}
\]

which is an exact infeasibility certificate for the endpoint row.

**[THEOREM]** The same one-edge class is obstructed on the whole incumbent
range, not only at the challenge point.  Positivity of the `atanh` series and
the certified rational `I3` upper endpoint give

\[
 \operatorname {atanh}(1/4)>
 \frac14+\frac{(1/4)^3}{3}=\frac{49}{192}
 >\frac{I_3}{2}.                                      \tag{3.6a}
\]

Consequently every `K<=I3/2` has `v=tanh K<1/4`, so the exact LP optimum
satisfies

\[
 p_*(K)=v^2<\frac1{16},\qquad
 \frac56-p_*(K)>\frac56-\frac1{16}=\frac{37}{48}.      \tag{3.6b}
\]

Thus `RC2-single-edge-conditional-domination` cannot improve or even reach the
incumbent anywhere in its entire relevant coupling range.

### 3.4 Independent finite-volume switching evaluation

**[COMPUTATION]** On the four-cycle `C4` embedded as a square plaquette of
`Z^3`, exact parity-subgraph enumeration at `v=6/25` gives

\[
 P_\varnothing=1+v^4=\frac{391921}{390625},
 \qquad P_{\{0,2\}}=2v^2=\frac{72}{625}.               \tag{3.7}
\]

**[COMPUTATION]** Direct enumeration of all `5^4` assignments of (3.1), with
both colour boundaries empty, gives

\[
 W_{\rm total}=P_\varnothing^2
 =\frac{153602070241}{152587890625},\qquad
 W_{0\leftrightarrow2}=P_{\{0,2\}}^2
 =\frac{5184}{390625}.                                 \tag{3.8}
\]

Therefore

\[
 \mathbb P(0\leftrightarrow2)
 =\frac{2025000000}{153602070241}
 =\left(\frac{P_{\{0,2\}}}{P_\varnothing}\right)^2,   \tag{3.9}
\]

verifying (3.2) exactly on the finite box.  The standalone verifier reconstructs
(3.8) by an edge-by-edge boundary/support dynamic program rather than the
producer's Cartesian-product enumeration.

**[UNRESOLVED]** The obstruction (3.6) says nothing about multi-edge block
events, a sharper self-contained percolation comparison, currents with sources,
or switching consequences not reducible to uniform single-edge iid domination.
Those are the remaining random-current directions after this finite class.

## 4. Classification and reproduction

**[COMPUTATION]** The artifact records eight universal-claim coverage rows:
`TARGET_BELOW_INCUMBENT`, `PM_PARSEVAL_BOUNDED`, `PM_SDP_OPTIMUM`,
`PM_POSITIVE_FLOOR_INFEASIBLE`, `RC_EDGE_WEIGHTS`, `RC_SWITCHING_C4`,
`RC_DOMINATION_LP`, and `PERC_PLANAR_5_6`.  The standalone verifier requires
this exact set and independently re-derives every decisive rational value.

**[COMPUTATION]** Run from the repository root:

```sh
.venv/bin/python experiments/e175_upper_endpoint4_paired.py
.venv/bin/python experiments/e176_upper_endpoint4_current.py
.venv/bin/python experiments/e177_upper_endpoint4.py
.venv/bin/python tests/test_upper_endpoint4.py
```

Each producer component and the combined producer prints a final `PASS` only
when all of its exact checks pass.  The combined producer writes
`results/bounds/upper_endpoint4.json`; the standalone verifier imports neither
producer component.
