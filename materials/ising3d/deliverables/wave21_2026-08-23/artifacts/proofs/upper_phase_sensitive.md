# Phase-sensitive Fourier convolution blocks strengthen MR4 but retain the zero uniform floor

**[COMPUTATION]** Artifacts: `experiment/e231_phase_sensitive_upper.py`,
`results/bounds/upper_phase_sensitive.json`, and
`tests/test_upper_phase_sensitive.py`.  The verifier imports no producer module.

## 0. Result and exact scope

**[THEOREM]** On the even cubic torus \(T_L=(\mathbb Z/L\mathbb Z)^3\), with
\(N=L^3\) and normalized amplitudes

\[
f_k=\frac1N\sum_{x\in T_L}\sigma_x e^{ik\cdot x},
\tag{0.1}
\]

the class \( \mathrm{PS4} \)-**convolution-block-Gram-L1** consists of all
\( \mathrm{MR4} \)-power-simplex-L1 rows plus, for every directed total momentum
\(k\), an Hermitian positive-semidefinite correlate of the block
\(x^{(k)}_q=f_qf_{k-q}\), subject to its prescribed power diagonal, involution
row equalities, and the sum-zero rows for \(k\ne0\).  This class is **strictly
stronger** than \( \mathrm{MR4} \)-power-simplex-L1: at \(L=4,K=6/25\), the
precisely specified point in Section 4 obeys every MR4 row but violates one
directed convolution polygon.

**[THEOREM]** Nevertheless the exact Green profile already used to close
MR4 lifts to this class for every even \(L\ge266\).  The all-size theorem
below gives exact polygon inequalities without a square-root or floating SDP
comparison.  Its zero mode satisfies

\[
0<\ p_0\le \frac{5168}{525L}
\tag{0.2}
\]

in the applicable coupling range, hence the new class has zero uniform
zero-mode floor in the thermodynamic limit.

**[UNRESOLVED]** The certified upper endpoint remains \(I_3/2\).  The result is
a split method theorem: one phase-sensitive row family genuinely rejects a
point that MR4 accepts, while the exact basis-case pseudo-moment is not
excluded.  Full cross-channel fourth-moment consistency, higher localizers,
DLR equations, and sourced or multi-edge current identities remain outside the
named class.

## 1. Exact convolution identity and derived four-point rows

**[LEMMA]** For every spin configuration and every \(k\),

\[
\sum_{q\in T_L^\ast} f_qf_{k-q}=\delta_{k,0}.
\tag{1.1}
\]

Indeed

\[
\sum_q f_qf_{k-q}
 =\frac1N\sum_x \sigma_x^2e^{ik\cdot x},
\]

because the factor \(1/N^2\) times exactly \(|T_L|=N\) ordered phase pairs
cancel for each \(x\), and \(\sigma_x^2=1\). \(\square\)

**[LEMMA]** At \(k=0\), equation (1.1) is Parseval's identity
\(\sum_kP_k=1\).  For \(k\ne0\), fix a phase pair
\(u,q\mapsto k-q\), sum the identity, and use \(P_qP_{k-q}\ge0\):

\[
\sqrt{P_0P_k}\le\sum_{q\notin\{0,k\}}\sqrt{P_qP_{k-q}},
\qquad k\ne0.
\tag{1.2}
\]

Averaging and using convexity gives

\[
2\sqrt{Q_{0,k}}
\le\sum_{q\notin\{0,k\}}\sqrt{Q_{q,k-q}},
\tag{1.3}
\]

where \(Q_{a,b}=\mathbb E[P_aP_b]\).  The configurationwise real polygon also
has exact power weights: on a regular involution orbit \(O=\{q,k-q\}\), both
representatives contribute \(\sqrt{P_qP_{k-q}}\), so the closed form is

\[
m_q\sqrt{P_qP_{k-q}}
\le\sum_{r\notin\{q,k-q\}}m_r\sqrt{P_rP_{k-r}},
\qquad m_r=\lvert\{r,k-r\}\rvert.
\tag{1.4}
\]

For every \(k\ne0\), this is a configuration-derived power-spectrum polygon
that no row of MR4 records.

The orbit for \(k=0\) is discarded from (1.3): it contributes the identity
\(\sqrt{Q_{0,0}}=\sqrt{Q_{0,0}}\), never an inequality on \(p_0\).

## 2. Directed convolution blocks

**[LEMMA]** Define

\[
x^{(k)}_q=f_qf_{k-q},\qquad
H^{(k)}_{q,r}=\mathbb E[x^{(k)}_q\overline{x^{(k)}_r}].
\tag{2.1}
\]

Every phase-sensitive fourth-spin contraction is an entry of one of these
blocks, because

\[
x^{(k)}_q\overline{x^{(k)}_r}
=\frac1{N^4}\sum_{x,y,u,v}
 \sigma_x\sigma_y\sigma_u\sigma_v
 e^{iq\cdot(x-y)+i(k-r)\cdot(u-v)}.
\tag{2.2}
\]

With the power matrix \(Q_{a,b}=\mathbb E[P_aP_b]\), the exact rows are

\[
\begin{aligned}
&H^{(k)}\succeq0,\qquad &H^{(k)}_{q,q}&=Q_{q,k-q},\\
&H^{(k)}_{k-q,r}&=H^{(k)}_{q,k-r},
\qquad &H^{(0)}&=Q,\quad Q\mathbf1=p,\\
&H^{(k)}\mathbf1&=0\quad(k\ne0),
\qquad &H^{(-k)}_{-q,-r}&=\overline{H^{(k)}_{q,r}},\\
&|H^{(k)}_{q,r}|^2&\le Q_{q,k-q}Q_{r,k-r}. &&
\end{aligned}
\tag{2.3}
\]

The PSD row is Gram positivity; row equality is
\(f_qf_{k-q}=f_{k-q}f_q\); the sum-zero row multiplies the configuration-wise
identity (1.1) by a second block entry.  The final row is Gram localizing.
No numerical optimization is involved.

**[THEOREM / equivalence lemma]** Fix \(k\ne0\), assume inversion symmetry
\(Q_{-a,b}=Q_{a,b}=Q_{a,-b}\), and set
\(O_q=\{q,k-q\}\), \(m_q=|O_q|\in\{1,2\}\), and
\(a_q=Q_{q,k-q}\ge0\).  The involution rows and sum-zero rows, together with
PSD, are equivalent to the complete power-polygon family

\[
m_q\sqrt{a_q}\ \le\
\sum_{r\notin\{q,k-q\}}m_r\sqrt{a_r}
\quad\text{for every }q.
\tag{2.4}
\]

*Proof.*  Take an involution-orbit representative \(q\), constants
\(\tau_q\in\{1,-1\}\), and set

\[
z_{k-q,s}=\tau_q z_{q,s}
\]

on a maximal set \(\{z_{q,s}:s=1,\dots,m_q\}\) of unit vectors with pairwise
dot products \(1-\tau_q\).  Define \(H\) semidefinately by

\[
H_{k-q,s;\,k-r,t}=H_{q,s;\,r,t}
=\tau_q\tau_r\sqrt{a_qa_r}\langle z_{q,s},z_{r,t}\rangle .
\tag{2.5}
\]

Here a noncentral orbit contributes the paired index sets for \(q\) and
\(k-q\), while a singleton orbit has one index.  The construction has the
prescribed diagonal and involution row.  Its row sum on the representative
orbit is

\[
\sum_{r,t}\tau_q\tau_r\sqrt{a_qa_r}\langle z,z'\rangle
=\tau_q\sqrt{a_q}\left(m_q\sqrt{a_q}
-\sum_{r\ne q}m_r\sqrt{a_r}\right).
\tag{2.6}
\]

Thus the sum-zero row vanishes exactly when (2.4) holds for that orbit.

Conversely, let \(H\) satisfy PSD, its prescribed diagonal, the involution
rows, and the sum-zero rows.  Because the left factor is \(m_q\), while the
configuration identity carries no factor of two, the configurationwise polygon
(1.4) requires \(m_r\), not the orbit count halved.  To recover it, define a
second block on involution-orbit representatives by

\[
J_{q,r}=\sqrt{a_qa_r}\langle z_{q,1},z_{r,1}\rangle .
\tag{2.7}
\]

Its exact sum-zero row is

\[
(J\mathbf1)_q
=m_q\sqrt{a_q}\left(m_q\sqrt{a_q}
 -\sum_{r\ne q}m_r\sqrt{a_r}\right).
\tag{2.8}
\]

Since \(J\) is a submatrix of twice the Gram block (2.5), its PSD is precisely
the polygon family (2.4); when equality fails there is no PSD realization of
the diagonal, involution rows, and sum-zero row simultaneously. \(\square\)

This theorem is the reason the named class, although written with PSD and
localizing rows, can be tested exactly on power supports without taking a
square root numerically.

## 3. Strict separation from MR4

**[THEOREM]** At \(L=4,N=64,K=6/25\), put \(\varepsilon=1/100\),

\[
p_0=\frac{97}{100},\qquad
p_{(2,0,0)}=p_{(0,2,0)}=p_{(0,0,2)}=\varepsilon,
\tag{3.1}
\]

and \(p_k=0\) otherwise; put \(Q=pp^T\).  Every MR4-power-simplex-L1 row is
satisfied, but the directed \(k=(2,0,0)\) convolution polygon is not.

**[LEMMA]** The point is feasible for the inherited power rows.

- Simplex: \(97/100+3/100=1\); inversion and cubic conjugacy are immediate.
- At an occupied nonzero mode \(\lambda=2\) and \(N\lambda=128\), so
  \[
  \frac1{2KN\lambda}=\frac{625}{768}
  =\frac{25/3}{192}>\frac1{100}.
  \tag{3.2}
  \]
- Every real-space correlation is
  \[
  p_0+\varepsilon\sum_{\text{three self-pair modes}}\cos(k\cdot z).
  \tag{3.3}
  \]
  One occupied mode contributes \(-1\) on an even coordinate and zero on the
  odd coordinates, hence every value lies in
  \(\{47/50,49/50,1\}\).
- At every unit edge direction, \(G(e)=49/50\).  The inherited finite energy
  row is
  \[
  G(e)\ge1-\frac{1-N^{-1}}{6K}
  =1-\frac{63/64}{36/25}
  =\frac{81}{256},
  \tag{3.4}
  \]
  which is strictly below \(49/50\).
- For \(Q=pp^T\), \(Q\mathbf1=p\), \(Q\ge0\), \(\operatorname{diag}Q\le p\),
  and the moment Schur covariance \(Q-pp^T\) is the zero matrix.

Thus every power-simplex level-one row passes exactly. \(\square\)

**[LEMMA]** At \(k=(2,0,0)\), the central orbit \(\{0,k\}\) is occupied and no
other ordered product \(p_qp_{k-q}\) is occupied: each \(q=k-r\) is the inverse
of \(r\) for this self-inverse \(k\), and the support is exactly
\(\{0,(2,0,0),(0,2,0),(0,0,2)\}\).  Therefore

\[
4Q_{0,k}=4p_0\varepsilon=\frac{97}{2500}>0,
\qquad
\sum_{r\notin\{0,k\}}\sqrt{Q_{r,k-r}}=0.
\tag{3.5}
\]

So the separating polygon is strict, not merely a feasibility deficit near
zero. \(\square\)

**[THEOREM]** Consequently
\(\mathrm{proj}_{p,Q}\mathrm{PS4}\subsetneq\mathrm{proj}_{p}\mathrm{MR4}\)
already at finite side \(L=4\).  The phase-sensitive class is not a relabeling
or repackaging of MR4. \(\square\)

## 4. Exact finite Green controls

**[COMPUTATION]** At \(L=4,6\), the producer recomputes

\[
p_0=1-\alpha C_L(0),\qquad
p_k=\frac{\alpha}{N\lambda(k)},\qquad
\alpha=\frac{25}{12},
\tag{4.1}
\]

using only rational cosine values
\(\cos(2\pi j/L)\).  It matches the inherited rational values

\[
p_0^{(4)}=\frac{1631}{9216},
\qquad
p_0^{(6)}=\frac{737509}{7185024}.
\tag{4.2}
\]

**[COMPUTATION]** Every nonzero directed momentum is covered.  There are
\(N-1\) momenta and, explicitly,

\[
\#\text{polygon rows}=
\frac{(N-1)N}{2}+4\left((L/2)^3-1\right).
\tag{4.3}
\]

For \(L=4\) this is \(2044\), and for \(L=6\) it is \(23776\).  For each
\(k\), the involution orbits partition all \(N\) indices, the central orbit is
present, and the producer hashes the per-momentum orbit counts.

**[THEOREM]** The stored Green powers satisfy every polygon row on
\(L=4,6\).

*Proof.*  Let \(m_q=|O_q|\le2\).  For every noncentral orbit,
\(\sqrt{a_q}<1/6\), because on the computed tori
\(\max_{k\ne0}p_k<1/36\).  There are at least \((N-4)/2\) nonnegative
remainder terms, each \((N-4)/2\) at least \(m_r\sqrt{a_r}\ge\sqrt{a_r}\) and
\(\sqrt{a_r}\ge p_{\min}\), with \(p_{\min}=\alpha/(6N)\).  Hence it is
sufficient that

\[
N\Bigl(\frac{N-4}{12}-p_{\max}\Bigr)
=\frac{N-4}{6}-\frac{2}{\lambda_{\min}}>0.
\tag{4.4}
\]

On \(L=4\), \(\lambda_{\min}=1\) and this margin is \(8\); on \(L=6\),
\(\lambda_{\min}=1/2\) and it is \(94/3\).

For the central orbit, the right side contains all \(N-2\) equal-orbit modes,
while the left squared by the polygon theorem is at most
\(4p_0\alpha/(N\lambda_{\min})\).  Therefore it is sufficient to compare exact
squares:

\[
\left(\frac{(N-2)\alpha}{6N}\right)^2
>
\frac{4p_0\alpha}{N\lambda_{\min}}.
\tag{4.5}
\]

For \(L=4\) the squared margin is
\(119575/1327104>0\); for \(L=6\) it is
\(19041025/172440576>0\).  Since every side of the unsquared polygon is
nonnegative, these exact squared inequalities prove the directed rows.
\(\square\)

## 5. All-size exact lift of the Green point

**[EXTERNAL]** Sections 4 and 5 of `proofs/upper_fourpoint.md` prove that, at
the benchmark-independent challenge \(K=6/25\), the Green branch is applicable
from even \(L\ge96\), and its exact zero mode rate is

\[
p_0=1-\frac{C_L(0)}{D_L}\le\frac{5168}{525L}.
\tag{5.1}
\]

At the endpoint \(K=I_3/2\), the same source uses
\(\alpha_L=\min(I_3^{-1},D_L^{-1})\).  For the uniform delayed bound used
below, replace \(I_3^{-1}\) by the larger \((2/L)\cdot9/8\), so that
\(\alpha=2/(3L)\) still obeys the finite infrared caps.  The certified
mpmath interval in `proofs/upper_infrared.md` has \(2(I_3/2)=I_3>16/27\),
hence \(2/(3L)<(8/9)I_3^{-1}\); equations (5.2)--(5.6) then give
\(\alpha\lambda(k)\ge8/(9L)\) at every nonzero mode.

**[LEMMA: dispersion lower bound]** For every even \(L\) and \(k\ne0\),

\[
\lambda(k)\ge\frac8{L^2}.
\tag{5.2}
\]

Indeed \(\lambda(k)=2\sum_i\sin^2(\pi k_i/L)\), and
\(\sin x\ge2x/\pi\) on \([0,\pi/2]\); at least one \(|k_i|_{\mathbb Z/L\mathbb Z}\)
is positive. \(\square\)

**[LEMMA: resistance diameter]** Write \(D_L=C_L(0)-\min_z C_L(z)\).  The
difference \(C_L(0)-C_L(z)\) is the effective resistance \(R_{\rm eff}(0,z)\)
on the unit-resistance cubic torus.  A path of at most \(3L/2\) unit edges
between the prefixed vertices has energy \(3L/2\), so

\[
D_L=\max_zR_{\rm eff}(0,z)\le\frac{3L}{2}.
\tag{5.3}
\]

*Proof.*  For the cycle coordinate alone, either pairwise-effective-resistance
physics or the explicit path voltage gives the bound.  Parallel edge
aggregates only reduce energy, and Thomson's principle sends any unit flow to
at most its energy.  Equality of \(C_L(0)-C_L(z)\) and resistance follows
because both solve the zero-sum Poisson equation. \(\square\)

**[LEMMA: profile branch]** For every \(0<K<I_3/2\), choose the Green point
with \(\alpha_L=1/(2K)\) on the finite base cutoff and \(\alpha_L=1/D_L\)
thereafter, exactly as in the inherited all-size theorem.  At \(K=I_3/2\), use

\[
\alpha_L=\min(I_3^{-1},D_L^{-1}).
\tag{5.4}
\]

In every case \(\alpha_L\le1/D_L\), so the GKS real-space rows hold and
\(\alpha_L\ge2/(3L)\); after the inherited rate cutoff, \(p_0\le5168/(525L)\).

**[THEOREM: noncentral polygons]** For every even \(L\ge4\) and every
noncentral orbit, first use \(\lambda(k)\le6\) to promote the infrared ceiling:

\[
\sqrt{a_q}
\le\sqrt{\frac{2}{3L}\,\frac{1}{2KN\lambda(k)}}
\le\frac1{3\sqrt{L}}<\frac16,
\tag{5.5a}
\]

for \(L\ge4\).  The remainder contains at least \((N-4)/2\) mode multiplicities
and each contributes \(m_r\sqrt{a_r}\ge\sqrt{a_r}\), with
\(a_r\ge(\alpha_L/(6N))^2\).  Thus

\[
\sum_{r\notin\{q,k-q,0,k\}}m_r\sqrt{a_r}
\ge\frac{N-4}{2}\,\frac{\alpha_L}{6N}.
\tag{5.5b}
\]

In normalized form, dividing \(m_q\sqrt{a_q}\) by its fullest possible
right-hand contribution \(m_q\alpha_L/(6N)\), a sufficient exact condition is

\[
\frac{N-4}{6}\ge\frac{2}{\lambda_{\min}}.
\tag{5.6}
\]

Using (5.2), the right side is at most \(L^2/4\), and

\[
\frac{L^3-4}{6}\ge\frac{L^2}{4}
\quad\Longleftrightarrow\quad
2L^3-3L^2-8\ge0,
\tag{5.7}
\]

which is strictly true for every \(L\ge4\). \(\square\)

**[THEOREM: central polygons]** For the orbit \(\{0,k\}\), the squared central
polygon follows whenever

\[
\frac{((N-2)\alpha_L/(6N))^2}
 {4p_0\alpha_L/(N\lambda_{\min})}
=
\frac{(N-2)^2\alpha_L\lambda_{\min}}{144p_0}>1.
\tag{5.8}
\]

Apply \(p_0\le5168/(525L)\), \(\alpha_L\ge2/(3L)\), and
\(\lambda_{\min}\ge8/L^2\).  It is sufficient that

\[
L\left(1-\frac2{L^3}\right)^2>27\frac{5168}{525},
\tag{5.9}
\]

or, defining

\[
F(L)=\frac{(L^3-2)^2}{L^5}-27\frac{5168}{525},
\tag{5.10}
\]

that \(F(L)>0\).  For real \(L\),

\[
F'(L)=1+\frac8{L^3}-\frac{20}{L^6}>0\qquad(L\ge2).
\tag{5.11}
\]

Directed rational arithmetic gives

\[
F(264)=-\frac{100019517465297}{56010902092800}<0,
\qquad
F(266)=\frac{1806844037121}{8323159178600}>0.
\tag{5.12}
\]

Therefore the first even side for this uniform central proof is \(266\), and
every even \(L\ge266\) follows by monotonicity. \(\square\)

**[THEOREM]** For every \(0<K\le I_3/2\), the Green pseudo-point is feasible
for \( \mathrm{PS4} \)-convolution-block-Gram-L1 for every even \(L\ge266\),
with the named-limitation witness extended by (5.1). ∎

## 6. Direction, finite floor, and limitation

**[LEMMA]** For \(p_k>0\), (1.3) rearranges to

\[
\sqrt{Q_{0,k}}
\le\sqrt{\mathbb E[P_0]\,\mathbb E[P_k]}=\sqrt{p_0p_k},
\qquad
p_0\le
\left(\frac{\sum_{q\notin\{0,k\}}\sqrt{Q_{q,k-q}}}{2\sqrt{p_k}}\right)^2
\tag{6.1}
\]

once \(Q=pp^T\), using Cauchy-Schwarz for the bridge.  The inequality is
therefore an **upper** bound on the zero-mode mass.  Its success on a Green
point is consistent with, and cannot itself refute, a vanishing uniform lower
floor.

**[LEMMA]** In this class an exact finite-\(L\), zero-\(p_0\) point cannot come
from its inherited pointwise GKS rows: summing \(G(z)\ge0\) gives

\[
p_0=\frac1N\sum_zG(z),
\tag{6.2}
\]

but does not force \(G(z)>0\).  The Euclidean minimization of
\(\sum_zG(z)\) subject to \(G(z)\ge0\) and \(G(0)=1\) has positive minimum only
if a further row is imposed.  As stored in the audited class with the exact
positive Green spectra and the finite computation, the relevant certified
point has \(p_0>0\).  More strongly, finite spin Parseval
gives physical \(M_L^2=\langle m^2\rangle\ge1/N\) for a uniform binary-spin
configuration measure, so \(M_L^2=0\) is impossible physically.  The
relaxation theorem below is about its uniform floor, not an exact finite zero.

**[THEOREM]** No positive \(\eta\), independent of \(L\), can be a class
certification of \(p_0=q_{00}\ge\eta\) for all sufficiently large even \(L\):
choose even \(L\ge266\) larger than \(5168/(525\eta)\), use Sections 4 and 5,
and note \(p_0\le5168/(525L)<\eta\). \(\square\)

**[THEOREM]** Because the polygon block rows admit the entire Green sequence,
the incumbent certified upper endpoint does not move:

\[
K_c\le\frac{I_3}{2}
\quad\text{remains the sharpened certified endpoint.}
\tag{6.3}
\]

The certified direction is unchanged: physical points lie in the relaxation,
so a positive class infimum would imply physical order; however the present
class infimum tends to zero.  The benchmark value \(K_c=0.221654626\) has no
selection, fitting, comparison, or validation role.

## 7. What is not claimed

**[UNRESOLVED]** The lift is a blockwise pseudo-moment construction.  It is
not asserted to come from a spin or Gibbs measure.  The class deliberately
omits the stronger data that might still exclude it:

- full equality of the complete fourth-moment tensor across all three
  representations of the same quartet and its full permutation symmetries;
- the diagonals and higher-degree localizers of
  \(|\hat\sigma|^4\) or higher Fourier powers;
- DLR conditional equations;
- sourced or multi-edge random-current identities.

Thus the proved zero-floor limitation is exact for the named class
\( \mathrm{PS4} \)-convolution-block-Gram-L1, not a statement about the
physical Ising magnetization and not a closure of every conceivable
phase-sensitive ancestor.

## 8. Reproduction

```text
nice -n 20 .venv/bin/python experiments/e231_phase_sensitive_upper.py
nice -n 20 .venv/bin/python tests/test_upper_phase_sensitive.py
```
