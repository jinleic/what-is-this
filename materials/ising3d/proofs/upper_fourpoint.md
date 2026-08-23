# An exact limitation theorem for a mode-resolved four-point infrared correction

**[COMPUTATION]** Artifacts: `experiments/e215_fourpoint_kernel.py`,
`experiments/e216_fourpoint_certificate.py`,
`experiments/e217_fourpoint_endpoint.py`,
`results/bounds/upper_fourpoint.json`, and
`tests/test_upper_fourpoint.py`.

## 0. Result and exact scope

**[THEOREM]** Let `MR4-power-simplex-L1` be the finite-torus certificate class
defined in Section 2.  It retains every Fourier-power mean and every entry of
their four-point second-moment matrix, together with the inherited infrared,
GKS, boundedness, Parseval, and energy rows.  At every finite side `L` and every
coupling `K`, its projection onto the two-point variables is *exactly* the
inherited two-point relaxation.  Consequently its exact optimum
magnetisation floor is exactly the two-point optimum.

**[THEOREM]** For every `0<K<=I3/2`, the exact optimum floors of this class
have a feasible sequence tending to zero as `L` tends to infinity through even
sides.  In particular, `MR4-power-simplex-L1` cannot certify a uniform positive
magnetisation floor and cannot improve the Watson endpoint.

**[COMPUTATION]** The benchmark-independent challenge is

\[
 K_\dagger=\frac6{25},\qquad 2K_\dagger=\frac{12}{25}.
 \tag{0.1}
\]

The certified rational lower endpoint `I3_lo` gives

\[
 I_3-2K_\dagger\ge I_{3,\mathrm{lo}}-\frac{12}{25}>0,
 \tag{0.2}
\]

so `K_dagger<I3/2` exactly.  The value `0.221654626` has no selection,
fitting, comparison, or validation role.

**[COMPUTATION]** On `T_12`, an exact rational primal point in the named class
has

\[
 p_0=M_L^2=
 \frac{94973361469787}{3762640205356032}<\frac1{32}
 \tag{0.3}
\]

with exact gap

\[
 \frac1{32}-p_0=
 \frac{22609144947589}{3762640205356032}>0.
 \tag{0.4}
\]

**[THEOREM]** The all-size obstruction is not inferred from (0.3).  For every
even `L>=316`, an analytically constructed exact rational torus kernel gives a
class-feasible point with

\[
 0\le p_0\le \frac{5168}{525L}<\frac1{32}.
 \tag{0.5}
\]

More generally the right side tends to zero, which rules out every fixed
positive floor, not only `1/32`.

**[UNRESOLVED]** The certified upper endpoint remains `I3/2`.  The theorem is
a limitation of the precisely named power-intensity level-one moment class.  It
does not exclude phase-sensitive Fourier contractions, higher localising
levels, a full spin moment matrix, DLR identities, or sourced/multi-current
switching inequalities.  Its rank-one counterpoints are pseudo-moments and are
not asserted to be physical Ising states.

## 1. Inherited walls and finite-torus normalisation

**[EXTERNAL]** `proofs/upper_infrared.md` proves the pointwise nonzero-mode
infrared ceiling

\[
 \widehat G_L(k)\le \frac1{2K\lambda(k)},\qquad
 \lambda(k)=3-\sum_{i=1}^3\cos k_i,
 \tag{1.1}
\]

and the endpoint `K_c<=I3/2`.  It also proves that the thermodynamic two-point
constraint system is saturated by the Green profile.  The factors of two and
the ordered-bond normalisation are inherited unchanged.

**[EXTERNAL]** `proofs/upper_endpoint4.md` proves exact obstructions for the
smaller aggregate class `PM4-aggregate-one-moment` and the one-edge current
class `RC2-single-edge-conditional-domination`; it explicitly leaves separate
mode variables open.  The class below addresses that open mode-resolved slot.

**[EXTERNAL]** `proofs/kc_upper_peierls.md` gives a valid but much weaker
Peierls endpoint and proves that its particular finite-head/geometric-tail
method cannot reach the incumbent.  No contour count or benchmark value is an
input here.

**[LEMMA]** On `T_L=(Z/LZ)^3`, with `N=L^3`, define the zero-mode-removed Green
kernel

\[
 C_L(z)=\frac1N\sum_{k\ne0}\frac{\cos(k\mathbin\cdot z)}{\lambda(k)}.
 \tag{1.2}
\]

It is the unique zero-sum solution of

\[
 \left(3-\frac12A\right)C_L(z)=\delta_{z,0}-\frac1N,
 \qquad \sum_zC_L(z)=0,
 \tag{1.3}
\]

where `A` is nearest-neighbour adjacency.

**[LEMMA]** Every value `C_L(z)` is rational.  Equation (1.3) is a linear
system over the rationals.  Its restriction to the zero-sum subspace is
invertible because its nonzero eigenvalues are the positive numbers
`lambda(k)`.  The unique solution therefore lies in `Q^N`.

**[LEMMA]** Put

\[
 m_L=\min_z C_L(z),\qquad D_L=C_L(0)-m_L.
 \tag{1.4}
\]

Then `C_L(0)>0`, `m_L<0`, and `C_L(z)<=C_L(0)`.  The last inequality follows
termwise from `cos(k.z)<=1` in (1.2).  The zero sum then forces a negative value,
so `D_L>C_L(0)>0`.

## 2. Derivation of the four-point inequalities before computation

**[LEMMA]** For a spin configuration, use the Fourier convention

\[
 \widehat\sigma(k)=\sum_{x\in T_L}\sigma_xe^{ik\cdot x},
 \qquad P_k(\sigma)=\frac{|\widehat\sigma(k)|^2}{N^2}.
 \tag{2.1}
\]

Configuration-wise Parseval and `sigma_x^2=1` give

\[
 P_k\ge0,\qquad P_{-k}=P_k,\qquad \sum_kP_k=1.
 \tag{2.2}
\]

**[LEMMA]** Define the mode-resolved two- and four-point variables

\[
 p_k=\mathbb E P_k=\frac{\widehat G_L(k)}N,
 \qquad Q_{k\ell}=\mathbb E[P_kP_\ell].
 \tag{2.3}
\]

The entries of `Q` are genuine four-spin quantities because

\[
 Q_{k\ell}=\frac1{N^4}\sum_{x,y,u,v}
 e^{ik\cdot(x-y)+i\ell\cdot(u-v)}
 \langle\sigma_x\sigma_y\sigma_u\sigma_v\rangle.
 \tag{2.4}
\]

**[LEMMA]** Equations (2.2)--(2.3) imply exactly

\[
 \begin{aligned}
 &p_k\ge0,\quad p_{-k}=p_k,\quad \sum_kp_k=1,\\
 &Q_{k\ell}\ge0,\quad Q\mathbf1=p,\quad
   \operatorname{diag}Q\le p,\\
 &\begin{pmatrix}1&p^T\\p&Q\end{pmatrix}\succeq0,
 \qquad Q_{-k,\ell}=Q_{k\ell}=Q_{k,-\ell}.
 \end{aligned}
 \tag{2.5}
\]

The row-sum identity follows by multiplying `P_k` by `sum_l P_l=1`; the
diagonal inequality follows from `0<=P_k<=1`; and the matrix is
`E[(1,P)^T(1,P)]`, hence positive semidefinite.  No numerical optimisation is
used to obtain these rows.

**[THEOREM]** The class `MR4-power-simplex-L1` consists of (2.5), the infrared
rows

\[
 p_k\le\frac1{2KN\lambda(k)}\quad(k\ne0),
 \tag{2.6}
\]

and the inherited real-space rows

\[
 0\le G_p(z):=\sum_kp_k\cos(k\cdot z)\le1,
 \quad
 G_p(e)\ge1-\frac{1-N^{-1}}{6K}.
 \tag{2.7}
\]

Every physical finite-torus Ising correlation supplies a feasible point, by
(1.1), (2.1)--(2.5), GKS-I, bounded spins, and the summed energy form of the
infrared ceiling.  Therefore a positive lower bound on the class objective
`p_0=M_L^2` is a valid physical lower bound; minimising over the relaxation has
the correct direction for an order certificate.

## 3. Exact projection and optimum theorem

**[THEOREM]** For every `p` satisfying the two-point rows (2.6)--(2.7) and the
first line of (2.5), set

\[
 Q=p\,p^T.
 \tag{3.1}
\]

Then `(p,Q)` satisfies every four-point row in (2.5).

**[LEMMA]** Elementwise nonnegativity is immediate.  Since `1^Tp=1`,
`Q1=p(p^T1)=p`.  Also `Q_kk=p_k^2<=p_k`.  Inversion symmetry is inherited
from `p`, and

\[
 \begin{pmatrix}1&p^T\\p&pp^T\end{pmatrix}
 =\binom1p\binom1p^{\!T}\succeq0,
 \qquad Q-pp^T=0.
 \tag{3.2}
\]

Thus (3.1) is an exact rational/algebraic rank-one certificate.

**[THEOREM]** If `B_L(K)` denotes the inherited two-point feasible set and
`M_L^{MR4}(K)` the class above, then

\[
 \operatorname{proj}_p M_L^{MR4}(K)=B_L(K),
 \qquad
 \inf_{(p,Q)\in M_L^{MR4}(K)}p_0
 =\inf_{p\in B_L(K)}p_0.
 \tag{3.3}
\]

One inclusion holds because the class definition retains all two-point rows;
the reverse inclusion is the lift (3.1).  Equation (3.3) is the exact optimum
statement.  It does not rely on a floating-point primal or dual SDP value.

**[THEOREM]** Any LP/SDP dual combination using only the named rows and claiming
`p_0>=eta` is exactly obstructed by any lifted two-point primal with
`p_0<eta`.  The four-point matrix cannot strengthen the dual because every
feasible two-point point has the zero-covariance lift (3.1).

## 4. Exact rational Green-ray counterpoints

**[LEMMA]** For fixed `K`, let

\[
 \alpha_L=\min\left\{\frac1{2K},\frac1{D_L}\right\},
 \quad p_0=1-\alpha_LC_L(0),
 \quad p_k=\frac{\alpha_L}{N\lambda(k)}\ (k\ne0).
 \tag{4.1}
\]

Then `p` satisfies every two-point row in the class.

**[LEMMA]** The probabilities are nonnegative and sum to one.  The cap (2.6)
holds because `alpha_L<=1/(2K)`.  Their inverse transform is

\[
 G_p(z)=1-\alpha_L\bigl(C_L(0)-C_L(z)\bigr).
 \tag{4.2}
\]

Since `0<=C_L(0)-C_L(z)<=D_L` and `alpha_L D_L<=1`, equation (4.2) lies in
`[0,1]`.  Summing the cap against `lambda` gives the energy row in (2.7).
Thus (4.1) is feasible, and (3.1) supplies its exact mode-resolved four-point
lift.

**[COMPUTATION]** The producer constructs `C_L` exactly for
`L=4,6,8,12`.  It stores signed-permutation orbit representatives and verifies
(1.3) at every site, zero spatial sum, rationality, the full inverse transform,
every cap, every row of (2.5), and the rank-one factorisation (3.2).  The finite
kernels are certificates rather than numerical approximations.

**[COMPUTATION]** At `L=12`,

\[
 \begin{aligned}
 C_{12}(0)&=\frac{733533368777249}{1567766752231680},\\
 m_{12}&=-\frac{39905507821}{3671584899840},\\
 D_{12}&=\frac{422619944041}{882751549680},\\
 \frac{12}{25}-D_{12}
 &=\frac{5503999027}{4413757748400}>0.
 \end{aligned}
 \tag{4.3}
\]

Hence the first branch of (4.1) applies with `alpha=25/12`; direct reduction
gives (0.3)--(0.4).

**[UNRESOLVED]** The `L=12` point alone has no thermodynamic consequence.  It
is only a finite exact countercertificate to the fixed row `p_0>=1/32` at that
side.  The all-size limitation comes from Section 5.

## 5. All-size rate and thermodynamic method limitation

**[EXTERNAL]** The exact lemmas in `proofs/mag_floor.md` give

\[
 C_L(0)\ge I_3-\frac{A}{L},\qquad
 A=\frac{992015}{408608},
 \tag{5.1}
\]

and

\[
 -m_L\le\frac{323}{200L},\qquad
 \frac{-m_L}{D_L}\le\frac{5168}{525L}.
 \tag{5.2}
\]

They are proved from a positive massive-walk kernel, exact lattice-sum bounds,
and rational upper bounds on `pi` and `sqrt(2)`; they are not empirical fits.

**[COMPUTATION]** At the rational challenge, directed rational arithmetic gives

\[
 I_{3,\mathrm{lo}}-\frac{A}{96}-\frac{12}{25}>0.
 \tag{5.3}
\]

Therefore every even `L>=96` has `D_L>=C_L(0)>=2K_dagger`, so the second branch
of (4.1) applies.

**[THEOREM]** For every even `L>=96`, the lifted point (4.1), (3.1) satisfies

\[
 p_0=1-\frac{C_L(0)}{D_L}
 =\frac{-m_L}{D_L}
 \le\frac{5168}{525L}.
 \tag{5.4}
\]

This is an all-size theorem based on exact inequalities, not a promotion of the
four computed sides.

**[COMPUTATION]** At `L=316`,

\[
 \frac1{32}-\frac{5168}{525\cdot316}
 =\frac{131}{1327200}>0.
 \tag{5.5}
\]

Thus (0.5) holds for every even `L>=316`.

**[THEOREM]** For any rational or real `eta>0`, choose an even

\[
 L\ge96,\qquad L>\frac{5168}{525\eta}.
 \tag{5.6}
\]

The exact class-feasible point then has `p_0<eta`.  Hence no `eta>0` and no
finite cutoff can make the class optimum uniformly at least `eta` at
`K_dagger`.

**[THEOREM]** The same conclusion holds for every `0<K<I3/2`: equation (5.1)
places the Green-diameter branch in force for all sufficiently large even `L`,
and (5.4) is independent of `K` thereafter.  At `K=I3/2`, the endpoint profile
`alpha_L=min(I3^{-1},D_L^{-1})` and the exact endpoint comparison in
`proofs/mag_floor.md` give the same `5168/(525L)` rate.  Combining this with
(3.3) proves the whole-incumbent-range method limitation stated in Section 0.

## 6. Direction and thermodynamic-limit passage

**[EXTERNAL]** For the ferromagnetic Ising model, if there are `eta>0` and
`L_0` such that the physical periodic-torus values obey
`M_L^2>=eta` for every even `L>=L_0`, then the plus state has positive
spontaneous magnetisation and `K_c<=K`.  This is the thermodynamic implication
used in the audited infrared proof.

**[LEMMA]** The direction can also be seen directly.  Since

\[
 M_L^2=\frac1N\sum_zG_L(z),
 \tag{6.1}
\]

a uniform positive average cannot be supported on any fixed finite set of
offsets as `N` tends to infinity.  For every fixed radius, some more distant
offset has a uniformly positive correlation along a subsequence.  Ferromagnetic
boundary monotonicity transfers periodic correlations to the plus state, and
the standard plus-state cluster limit identifies the limiting correlation with
`m_+(K)^2`.

**[THEOREM]** Because every physical point belongs to the relaxation, a class
lower bound `inf p_0>=eta>0` would imply the physical premise above and hence a
valid strict upper endpoint.  Equations (3.3) and (5.4) show that the named
class cannot produce that premise.  The inequality direction is therefore:
class infimum `<=` physical `M_L^2`; positivity of the former would certify
order, while a small pseudo-moment only proves a method limitation.

## 7. Counterexample scope

**[THEOREM]** The counterpoint uses a Dirac measure on the Fourier-power
simplex: `Q=pp^T`.  This is exactly why every level-one power moment row is
saturated.  It is a legal point of the named relaxation but need not arise from
any spin configuration or Gibbs measure.

**[UNRESOLVED]** The class discards Fourier phases.  In particular it does not
impose the configuration-wise convolution identities

\[
 \sum_k\widehat\sigma(k)\widehat\sigma(q-k)
 =N^2\delta_{q,0},
 \tag{7.1}
\]

which encode all local equations `sigma_x^2=1`, rather than only their global
Parseval sum.  Phase-sensitive contractions can exclude the rank-one power
pseudo-measure and are not covered by the theorem.

**[UNRESOLVED]** Full spin moment matrices, higher simplex localisers, DLR
conditional equations, sourced random-current identities, and multi-edge
switching events also lie outside `MR4-power-simplex-L1`.  No negative claim is
made about those stronger classes.

**[UNRESOLVED]** No assertion is made that the true physical `M_L^2` tends to
zero at the challenge.  Only the optimum of the explicitly listed relaxation
is shown to have vanishing floor.

## 8. Reproduction and independent verification

**[COMPUTATION]** Run the producers from the repository root:

```text
nice -n 10 .venv/bin/python experiments/e215_fourpoint_kernel.py
nice -n 10 .venv/bin/python experiments/e216_fourpoint_certificate.py
nice -n 10 .venv/bin/python experiments/e217_fourpoint_endpoint.py
```

**[COMPUTATION]** The combined producer writes a top-level
`meta / data / checks` artifact.  Every stored passing check is evaluated from
exact arithmetic.  The standalone verifier imports no producer module and
rechecks the rational Poisson systems, orbit partitions, exact algebraic mode
sums, moment rows, challenge gaps, and all-size rational cutoffs.

**[COMPUTATION]** The independent verification command is:

```text
nice -n 10 .venv/bin/python tests/test_upper_fourpoint.py
```
