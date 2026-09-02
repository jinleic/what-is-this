# Liu H2 channel theory: a dimension-free rank-one decomposition

## Outcome

The useful result is an exact functional decomposition, not a six-dimensional
cover.  Put

\[
 \phi(s)=\sqrt{Q_2(s,s)},\qquad
 A(s,t)=P_2(s,t)+\phi(s)\phi(t),\qquad
 B(s,t)=Q_2(s,t)-\phi(s)\phi(t).
\]

Then `R=A+B`, and exact expansion gives

\[
\boxed{
 T(\mu,\nu)=\langle\mu,B\mu\rangle+\langle\nu,B\nu\rangle
 +2\langle\mu,A\nu\rangle
 +(\langle\mu,\phi\rangle-\langle\nu,\phi\rangle)^2.}
\tag{1}
\]

Consequently the scalar sandwich

\[
 0\le B(s,t)=Q_2(s,t)-\sqrt{Q_2(s,s)Q_2(t,t)}\le R(s,t)
 \quad(0\le s,t\le1)
\tag{2}
\]

proves Lemma C for **all** finite equal-mass nonnegative measures.  The left
half of (2) is `B>=0`; the right half is `A=R-B>=0`.  This is the only
remaining analytic/certification input.  Its two Arb certificates are owned by
`liu9_h2_copositive.py` and `liu9_h2_twovar_lemmas.py`, respectively; this
file intentionally does not duplicate their box work.

| route | verdict | decisive number or formula |
|---|---|---|
| 1. Conditional positive definiteness of `K` | **REFUTED** | `d=delta_(3/8)-delta_(7/8)` gives certified `beta<d,Kd> = -0.027550630591841221502363088586120520039...` |
| 2. Diananda `N+S` | **PROVED-CONDITIONAL-ON-THE-SANDWICH** | `S=vv^T`, `N_ij=A(x_i,y_j)+A(x_j,y_i)+B(x_i,x_j)+B(y_i,y_j)` |
| 3. Bernstein/Mercer truncation | **OPEN-WITH-OBSTRUCTION** | rigorous additive tail is `epsilon_T(n)` in (10); `epsilon_T(4096)<=0.326543475473237`, while `inf T=0` |
| 4. Pointwise `kappa` | **OPEN-WITH-OBSTRUCTION** | sharp candidate `kappa_0=4m beta/(2m-1)=1.05313420857430436016...<2`; the global upper bound and its measure-level tensorization are not proved here |

The exact symbolic identities are implemented at
`liu9_h2_channel_analysis.py:99-129`; the measure and matrix forms are emitted
at `liu9_h2_channel_analysis.py:549-554`.

## 1. The rank-one Diananda decomposition

For paired supports `(x_i,y_i)`, define

\[
 v_i=\phi(x_i)-\phi(y_i)
\]

and

\[
 N_{ij}=A(x_i,y_j)+A(x_j,y_i)+B(x_i,x_j)+B(y_i,y_j).
\]

Substitution of `P2=A-phi tensor phi` and `Q2=B+phi tensor phi` gives, with
zero residual in free symbols,

\[
 (M_T)_{ij}=N_{ij}+v_i v_j. \tag{3}
\]

If (2) holds, `N` is entrywise nonnegative and `S=vv^T` is PSD of rank at
most one.  This is stronger than the order-three consequence of Diananda:
(3) works in every order and (1) works for arbitrary finite measures.
The deterministic discovery check used 4,096 three-atom pairs and obtained

* maximum absolute identity residual `5.5511151231257827e-17`;
* minimum sampled `A = 2.6904219846896049e-6`;
* minimum sampled `B = -1.3877787807814457e-17` (roundoff on its exact zero);
* minimum sampled `N_ij = 1.1497723889572486e-5`.

These float64 values are evidence only.  The algebraic residual in (3) is
exact.

### 1.1 Why the less careful decompositions fail

The obvious split

\[
 N^{(R)}_{ij}=R(x_i,y_j)+R(x_j,y_i),\qquad
 S^{(K)}_{ij}=\beta[K(x_i,x_j)+K(y_i,y_j)-K(x_i,y_j)-K(x_j,y_i)]
\]

has entrywise `N^(R)>=0`, but `S^(K)` is not PSD.  At

\[
 (x_i,y_i)=((0,1),(0,1),(1/2,0))
\]

its smallest eigenvalue is `-0.07181059493033233`.

Nor may every positive entry simply be moved to `N`.  At

\[
 (x_i,y_i)=((x_\star,0),(x_\star,0),(0,x_\star))
\]

write `d=beta K(xstar,xstar)>0`.  The sign-clipped PSD candidate has exact
smallest eigenvalue

\[
 d(1-\sqrt2)=-0.028683036977071422\ldots. \tag{4}
\]

The positive edge discarded by clipping is precisely what completes the
rank-one square.  Formula (3), rather than sign clipping, keeps that edge in
`S`.

### 1.2 Is `phi=sqrt(Q2 diagonal)` forced?

For any rank-one ansatz with both pointwise pieces nonnegative, its diagonal
must satisfy

\[
 \max\{0,-P_2(s,s)\}\le \phi(s)^2\le Q_2(s,s). \tag{5}
\]

At `s=xstar`, `R(xstar,xstar)=0`, so `-P2=Q2` and (5) forces
`phi(xstar)^2=Q2(xstar,xstar)`.  It also forces `phi(0)=phi(1)=0`.
Thus **no choice of phi can remove the tight interior zero**.

The diagonal constraints do not force it away from these points.  Whether all
off-diagonal constraints restore uniqueness is a separate question.  If
canonical `phi` is changed only at `s=1/2` by a scale `r`, deterministic
discovery optimization of all partners `t` gives

\[
 r_{\rm offdiag}\ge0.918282886810733,
 \qquad r_{\rm diag}\ge0.943952927473383.
\]

Hence `r=0.95` retains both pointwise inequalities in the discovery check;
the off-diagonal maximum occurs at `t=0.599610592485236`.  The canonical
value is `phi(1/2)=0.24928153094486544`.  Thus nonuniqueness away from
`xstar` is supported numerically, not certified here.  The rigorous design
point is nevertheless fixed by (5): even a better admissible `phi` cannot
turn (2) into a strict global inequality because (5) is an equality at
`xstar`.

## 2. Route 1 — `K` is not conditionally positive definite

### Verdict: REFUTED

Take the exact signed measure

\[
 d=\delta_{3/8}-\delta_{7/8},\qquad d([0,1])=0.
\]

At both 320 and 512 Arb bits,

\[
 \langle d,Kd\rangle
 =-0.2753615762514741428518548157942339566375095606\ldots<0,
\]
\[
 \beta\langle d,Kd\rangle
 \in[-0.0275506305918412215023630885861205200390352255\ldots].
\]

with ball radius `4.98e-69`; specifically, the Arb **upper endpoint is
negative**.  This is the correct direction: the exact value is at most its
upper endpoint, so `upper<0` proves negativity.  The two mandatory mutations
also fail at both precisions:

* flipping the channel sign has upper endpoint `+0.0275506305918412...`;
* deleting the `-2K(s,t)` cross term has upper endpoint
  `+0.102466283119593...`.

The Arb code and the directional guards are at
`liu9_h2_channel_analysis.py:132-208`.

The certified rational witness is not the continuum minimum.  A 100-dps
boundary solve gives the discovery candidate

\[
 (s,t)=(0.33947854122264469790554106822801418744\ldots,1),
\]
\[
 \langle d,Kd\rangle=-0.83265720781258145276283391850517167731\ldots,
 \qquad c=-0.08330948512993773158095257336665329083\ldots. \tag{6}
\]

### 2.1 Fine-grid centred spectrum

On the 801-point uniform grid, with `H=I-11^T/801`, float64 discovery gives

\[
 \lambda_{\min}(HKH)=-40.748115129507809,
 \qquad \lambda_{\min}/801=-0.050871554468798763. \tag{7}
\]

There are 20 eigenvalues below `-1e-10`, two above `1e-10`, and the bottom
eigenvector has two sign changes.  After normalizing its positive and negative
parts to unit mass:

* positive part: mean support `0.886013944390027`, effective support size
  `175.4847`, largest weight `0.00957872077351324` at `1`;
* negative part: mean support `0.360453571264291`, effective support size
  `436.9555`, largest weight `0.00283408043429008` at `0.3525`;
* resulting `K` energy `-0.32545788645727491`, beta-scaled channel
  `-0.032562894667617172`.

### 2.2 The proposed midpoint-convexity characterization is incomplete

Put `w_ab=log(1-pi(a,b))` and `g(w)=h(1-exp(w))`.  Then

\[
 K(s,s)+K(t,t)-2K(s,t)
 =[g(w_{ss})+g(w_{tt})-2g(\bar w)]
 +2[g(\bar w)-g(w_{st})],\quad
 \bar w=(w_{ss}+w_{tt})/2. \tag{8}
\]

The second term generally does not vanish.  Exact algebra gives

\[
 [1-\pi(s,t)]^2-[1-\pi(s,s)][1-\pi(t,t)]
 =(s-t)^2\big((s+t-1)^2+1-s^2t^2\big)\ge0,
\]
so `w_st>=bar_w`, with strict inequality off the diagonal.  At `(3/8,7/8)`:

* `w_st-bar_w = 0.4238300301288183113...`;
* midpoint term `-0.3383510509044814641...`;
* displacement term `+0.06298947465300732129...`;
* total `-0.27536157625147414285...`.

More decisively, for the boundary minimizer `t=1`, `bar_w=-infinity`; the
midpoint term is positive while the displacement term makes the total
negative.  Thus failure is not “precisely where `g` fails convexity.”  The
curvature of `g` changes at probability
`0.78218829428019990122029707592674478018...`, or
`w=-1.52412432465752932602616841324896354586...`, but (8) has an independent
cross-log displacement.

## 3. Route 2 — exact `N+S`

### Verdict: PROVED-CONDITIONAL-ON-THE-SANDWICH

Equations (1)--(3) are the requested explicit, support-only natural
Diananda decomposition.  They avoid six-dimensional boxing and are
uniform in the number of atoms.  Route 2 closes Lemma C as soon as the two
scalar halves of (2) are certified.  Until those separately owned artifacts
land, the honest status is conditional rather than an unconditional theorem.

The difficult half `B>=0` is equivalently

\[
 h(\pi(s,t))^2\ge h(\pi(s,s))h(\pi(t,t)). \tag{9}
\]

The protocol arguments go in the *opposite* Cauchy--Schwarz direction:

\[
 \pi(s,s)\pi(t,t)-\pi(s,t)^2=s^2t^2(s-t)^2\ge0.
\]

Therefore (9) is genuinely an entropy-shape statement; monotonicity of `h`
or ordinary Cauchy--Schwarz cannot prove it.  The 801-grid check found only
roundoff `-1.39e-17` on exact zero strata.  The other half `A>=0` has the
forced nondegenerate zero at `(xstar,xstar)`.

## 4. Route 3 — Bernstein/Mercer truncation

### Verdict: OPEN-WITH-OBSTRUCTION

Use the tensor Bernstein operator `B_n` with `n>=4`.  Since
`|partial_s pi|,|partial_t pi|<=1`,

\[
 \mathbb E|X/n-s|\le\frac1{2\sqrt n},\qquad
 \mathbb E|Y/n-t|\le\frac1{2\sqrt n}.
\]

The binary-entropy continuity inequality
`|h(u)-h(v)|<=h(|u-v|)` and Jensen give the rigorous uniform tails

\[
 \epsilon_K(n)=h(1/\sqrt n),
\]
\[
 \epsilon_P(n)=(1-\beta)h(1/\sqrt n)
 +\frac{h(1/(2\sqrt n))}{m}
 +\frac{\log2}{2m\sqrt n},
\]
\[
 \epsilon_R(n)=h(1/\sqrt n)
 +\frac{h(1/(2\sqrt n))}{m}
 +\frac{\log2}{2m\sqrt n},
\]
\[
 \boxed{\epsilon_T(n)=2h(1/\sqrt n)
 +\frac{2h(1/(2\sqrt n))}{m}
 +\frac{\log2}{m\sqrt n}.} \tag{10}
\]

For each displayed Arb bound, the direction is

\[
 |T-B_nT|\le\text{exact expression}\le
 \operatorname{upper}(\text{Arb expression}).
\]

At 320 bits the certified upper expressions are:

| `n` | upper tail expression |
|---:|---:|
| 64 | `1.6513774753623700121322487553542...` |
| 256 | `0.9883155838425295602569586375105...` |
| 1024 | `0.5739797084646206472666988214567...` |
| 4096 | `0.3265434754732365143987317876055...` |

The first orders at which (10) falls below fixed tolerances are

| tolerance | first `n` |
|---:|---:|
| `1e-1` | 67,540 |
| `1e-2` | 12,569,023 |
| `1e-3` | 1,991,043,446 |
| `1e-6` | 5,087,256,357,468,862 |

No finite order suffices globally with this additive tail: `T` has exact
zeros, so its global positive margin is zero.  A finite truncation plus any
positive additive remainder cannot prove nonnegativity at those zeros without
exact local factorization.  The straightforward PSD relaxation also fails on
the 801-grid Bernstein coefficient matrix:

* `lambda_min(R_grid) = -12.468076358775775`;
* `lambda_min((2Q-R)_grid) = -0.012949302032404495`.

Thus the truncation still leaves a high-order copositivity problem and a
zero-local remainder problem; it is not an analytic shortcut.

## 5. Route 4 — sharp pointwise comparison

### Verdict: OPEN-WITH-OBSTRUCTION

Let

\[
 c_2(s,t)=\beta[K(s,s)+K(t,t)-2K(s,t)].
\]

On the boundary `t=1`, exact cancellation gives

\[
 R(s,1)=a h(s),\qquad a=1-\frac1{2m},
\]
\[
 -\frac{c_2(s,1)}{R(s,1)}
 =\frac{2\beta}{a}-\frac{\beta K(s,s)}{a h(s)}.
\]

Because `K(s,s)/h(s)->0` as `s->0+`, every valid pointwise constant must
satisfy

\[
 \boxed{\kappa\ge\kappa_0:=\frac{2\beta}{a}
 =\frac{4m\beta}{2m-1}
 =1.05313420857430436016386400743722888219\ldots.} \tag{11}
\]

Arb at both 320 and 512 bits encloses (11) in

`[1.05313420857430436016386400743722888219194036 +/- 3.33e-45]`,

so the **upper endpoint is below 2**.  The direction is important: this proves
`kappa_0<2`, but it does not prove that no interior pair has a larger ratio.
An 801-grid discovery maximum is `1.0507939414519576`; the singular boundary
limit is missed by every fixed grid.  Along `s=10^-16,t=1`, the deficit from
(11) is `4.88277437717519294e-15`.

The interior binding point is not sharp.  Arb gives

\[
 K_{12}(x_\star,x_\star)
 =-1.59260964938831675827994499403826726710\ldots
\]

and the exact quadratic ratio

\[
 \frac{-4\beta K_{12}}{R_{11}-R_{12}}
 =0.51007705032962250533090955868792334334\ldots, \tag{12}
\]

strictly below (11).  Thus the earlier expectation that `xstar` controls
`kappa` is refuted numerically and by the certified comparison of (11) and
(12); the singular corner `(0,1)` controls the observed supremum.

### 5.1 Why the two-point bound does not automatically tensorize

In the rank-one notation,

\[
 c=\langle\mu,B\mu\rangle+\langle\nu,B\nu\rangle
 -2\langle\mu,B\nu\rangle+(\Phi_\mu-\Phi_\nu)^2,
 \quad \Phi_\mu=\langle\mu,\phi\rangle.
\]

If the point-mass comparison is integrated against `mu tensor nu`, the
difference between the desired general expression and that integrated
pointwise expression is

\[
 [\langle\mu,B\mu\rangle-\operatorname{Var}_\mu(\phi)]
 +[\langle\nu,B\nu\rangle-\operatorname{Var}_\nu(\phi)]. \tag{13}
\]

It has no fixed nonnegative sign.  For
`lambda=(delta_0+delta_(1/2))/2`, `B(0,1/2)=0` and

\[
 \langle\lambda,B\lambda\rangle-
 \operatorname{Var}_\lambda(\phi)
 =-0.01553532041755397813293779310888594327\ldots. \tag{14}
\]

So naive integration is genuinely invalid.  This is an obstruction to the
proof route, not a counterexample to the general inequality.  A deterministic
4,096-draw, three-atom discovery search found maximum general ratio
`0.72305907255120427`, while the explicit family

\[
 \mu=(1-p)\delta_0+p\delta_s,\qquad \nu=\delta_1
\]

has the exact ratio

\[
 \kappa_0-\frac{p\beta K(s,s)}{a h(s)}
\]

and approaches `kappa_0` from below as `p->0`.  No example exceeding
`kappa_0` was found.  The statement that the same constant controls all
measures therefore remains **open**, for the precise tensorization reason
(13)--(14).

## 6. Reproduction and soundness checks

Run from the project root:

```bash
nice -n 19 ./.venv/bin/python -B uc/liu9_h2_channel_analysis.py > /tmp/channel-a.json
nice -n 19 ./.venv/bin/python -B uc/liu9_h2_channel_analysis.py > /tmp/channel-b.json
cmp /tmp/channel-a.json /tmp/channel-b.json
shasum -a 256 /tmp/channel-a.json
```

The commands actually run for this revision produced:

* `cmp`: exit 0;
* stdout file SHA-256:
  `6e3b02ec8b050a6aff4f300c3162a033da4f20061c8b5849d01e40d1da00b2e5`;
* byte count: `25,212`;
* payload-internal SHA-256 (hash before adding its own field):
  `5e5d62bf7504fd49bd6f7e935d02f3b9fdd2943ab5d9c1d9247a1bb93b8d51fa`.

The report contains four mandatory mutant rows (two mutations at each of 320
and 512 bits), and all four fire.  No wall-clock time, randomness without a
fixed seed, or machine-dependent timestamp enters the JSON payload.
