# CCF self-similar profile — PAUSED technique validation

> **PAUSED, and not on the canonical path.** The canonical target is Clay Navier–Stokes,
> Fefferman (C)/(D): see [`../ns/README.md`](../ns/README.md).
>
> This subproject is now known to be **structurally** off that path, not merely
> unbridged. Lemma 1 (`../ns/README.md`) forces any exact self-similar NS blowup to have
> $\beta=1/2$, and Tsai's Theorem 2 (ARMA 143 (1998) 29–51) closes $\beta=1/2$ assuming
> only distributional NS plus a local energy bound. So the entire profile-based program
> — this, Chen–Hou, and the DeepMind unstable-singularity line — cannot reach the Clay
> statement by construction, whatever its merits as inviscid model mathematics.
>
> Retained for two things it produced: the exactly-invertible
> [`tail_basis.py`](tail_basis.py) family, and the lesson recorded below about
> validating an operator on the decay class you actually need.

**Prize relevance: none.** 1D, inviscid, whole-line. VISCOSITY: no — the candidates are
inviscid or fractional $\alpha\le1$; constant $\nu\Delta$ ($\alpha=2$) lies far outside.
DOMAIN: no — 1D $\mathbb{R}$, not $\mathbb{R}^3$/$\mathbb{T}^3$. See
[`../docs/GAP_MAP.md`](../docs/GAP_MAP.md).

## The equation

Córdoba–Córdoba–Fontelos, J. Math. Pures Appl. 86 (2006) 529–540, eq (1):

$$f_t - (Hf)\,f_x = 0, \qquad Hf(x) = \frac{1}{\pi}\,\mathrm{PV}\!\int_{\mathbb{R}}\frac{f(y)}{x-y}\,dy, \qquad \widehat{Hf}(\xi) = -i\,\mathrm{sgn}(\xi)\,\hat f(\xi).$$

Vorticity form ($\omega=f_x$, $u=Hf$, $u_x=H\omega$), arXiv:2511.22819 eq (2):

$$\omega_t - u\,\omega_x = \omega\,u_x, \qquad u(0)=0.$$

Self-similar ansatz, eq (3): $\omega=(1-t)^{-1}\Omega(y)$, $y=x(1-t)^{-(1+\lambda)}$,
$U(y)=\int_0^y H\Omega$.

## Profile equation — derived here, cross-checked

Substituting gives $\omega_t=(1-t)^{-2}[\Omega+(1+\lambda)y\Omega']$,
$\omega_x=(1-t)^{-(2+\lambda)}\Omega'$, $u=(1-t)^{\lambda}U$, hence
$u\omega_x=(1-t)^{-2}U\Omega'$ and $\omega u_x=(1-t)^{-2}\Omega\,H\Omega$, so

$$\textbf{(P)}\qquad \Omega + \big[(1+\lambda)y - U\big]\Omega' - \Omega\,(H\Omega) = 0.$$

This reproduces arXiv:2511.22819 eq (4) — an independent check on all sign and
scaling conventions before any numerics were written.

### Exact first integral

Put $G=(1+\lambda)y-U$, so $G'=(1+\lambda)-H\Omega$. Then

$$(G\Omega)' = G'\Omega + G\Omega' = (1+\lambda)\Omega - \Omega H\Omega + \big(-\Omega + \Omega H\Omega\big) = \lambda\,\Omega,$$

using (P) for $G\Omega'$. With $F(y)=\int_0^y\Omega$ and $F(0)=0$, and noting
$G(0)\Omega(0)=0$:

$$\textbf{(Q)}\qquad \big[(1+\lambda)y - U(y)\big]\,\Omega(y) = \lambda\,F(y).$$

**(Q) is algebraic in $\Omega$** — no derivative of the unknown appears. That is what is
discretized: it removes the worst source of conditioning from the collocation system.
Elementary, and plausibly known, but not stated in any source retrieved for this project.

A further consequence, useful later: $U = HF$ exactly, since $U'=H\Omega=H(F')=(HF)'$ and
both vanish at $0$ ($F$ even ⇒ $HF$ odd).

## Structure

- **Parity.** (P) is consistent only for $\Omega$ **odd**; then $H\Omega$ and $F$ are even,
  $U$ and $G$ odd. (Consistent with $f$ even and positive, $\omega=f_x$ odd, $u(0)=0$.)
- **Symmetry.** The only continuous symmetry is the dilation $\Omega(y)\mapsto\Omega(by)$.
  An amplitude rescale $\Omega\mapsto a\Omega$ is **not** a symmetry: substituting
  $a\Omega(by)$ into (P) yields the original equation multiplied by $a$ only if $a=1$.
  Normalizing $\Omega'(0)=1$ therefore removes the last degree of freedom and makes
  $\lambda$ an isolated nonlinear eigenvalue.
- **Far field.** Balancing $\Omega+(1+\lambda)y\Omega'\approx0$ gives
  $\Omega\sim C\,\mathrm{sgn}(y)|y|^{-p}$ with $p=1/(1+\lambda)$. For $\lambda\approx1.18$,
  $p\approx0.4586$ — very slow decay, and the source of every numerical difficulty here.
  Accordingly we solve for the bounded $V$ in
  $\Omega(y)=y(1+y^2)^{-(p+1)/2}V(y)$, $V$ even, $V(0)=1$, $V(\infty)=C$.

## Target constants

From arXiv:2509.14185 Fig. 2 (CCF column) `[REPORTED]`:

| mode | $\lambda$ | $\log_{10}$ max residual (Fig. 3) |
|---|---|---|
| stable | **1.1807776628998** | −13.714 |
| 1st unstable | 0.6057337012032 | −13.589 |
| 2nd unstable | 0.4713248638620 | −6.664 |

Corroboration: Eggers–Fontelos (Nonlinearity 2019) give stable $\lambda\approx1.18078$.

**Recorded conflict.** Two sweep agents reported the 2nd unstable value as
$\lambda_2=0.4703$ with residual $O(10^{-7})$ from the paper text, against
$0.4713248638620$ from Fig. 2. Unresolved; do not rely on either until the paper is
re-read directly. Not needed for the current stable-mode target.

Related tables retained for later: IPM $\lambda$ = 1.0285722760222 (stable),
0.4721297362414, 0.3149620267088, 0.2415604743989; Boussinesq = 1.9205599746927,
1.3990961221852, 1.2523481636489, 1.1842500861997, 1.1448833556188 (last unvalidated).
The Boussinesq list is corroborated by the paper's own empirical fit
$\lambda_n\approx 1/(1.4187n+1.0863)+1$, which reproduces all five to ~3 decimals; a
sweep agent that omitted the 1.2523 row is therefore the one in error.

## Method

1. Compactify with $y=\tan(\theta/2)$ on a **half-cell-offset** uniform $\theta$ grid,
   $\theta_j=-\pi+2\pi(j+\tfrac12)/n$. Offsetting is mandatory: a node at $\theta=-\pi$
   makes $\tan(\theta/2)\sim10^{16}$ and destroys the transform by cancellation.
   The grid is antisymmetric-paired, so parity constraints are exact.
2. Hilbert transform via the Cayley map (`hilbert.py`):
   $Hf(\tan(\theta/2)) = H_{\text{circ}}[g](\theta) + C[f]$, with
   $C[f]=-\tfrac1\pi\int_{\mathbb{R}} f(x)\,x/(1+x^2)\,dx$, evaluated by FFT.
3. Discretize **(Q)** on the positive half-grid; unknowns are $V$ plus $\lambda$;
   normalization $V(0)=1$; solve by damped least squares.
4. Validate every component against exact or `mpmath` high-precision references before
   trusting any $\lambda$.

## Status

Working, with one component defect under repair.

- `hilbert.py` — Hilbert operator. Validated to **machine precision** on three exact
  transform pairs and on $H^2=-I$ (8.96e-16) for **rapidly decaying** data.
- **Defect found (this is the current work).** For **algebraically decaying** data the
  additive constant $C[f]$ is computed by a uniform-grid quadrature of an integrand
  singular like $(\pi-\theta)^{p-1}$, giving only $O(N^{-(1-p)})$ accuracy —
  measured $O(N^{-1/2})$: errors 3.34e-2 → 1.77e-2 → 9.38e-3 at $N=1024,4096,16384$,
  and nearly independent of $y$, i.e. a pure constant offset.
  This is fatal rather than cosmetic: a constant error $\delta$ in $H\Omega$ shifts $U$
  by $\delta y$, which is **exactly equivalent to shifting $\lambda$ by $\delta$**.
  The first solve consequently ran $\lambda$ to its guard rail.
  The earlier validation missed it because Poisson-kernel test data decays like $y^{-2}$
  — too fast to excite the singular endpoint behaviour.
- **Fix in progress.** Split off the model tail $m(y)=y(1+y^2)^{-(p+1)/2}$, whose
  constant is exact:
  $C[m]=-\tfrac1\pi\int_{\mathbb{R}} y^2(1+y^2)^{-(p+3)/2}dy = -\tfrac1\pi B(\tfrac32,\tfrac p2)$,
  from $\int_0^\infty x^{a-1}(1+x^2)^{-b}dx=\tfrac12 B(\tfrac a2, b-\tfrac a2)$ with
  $a=3$, $b=(p+3)/2$. Quadrature then acts only on the fast remainder.
- Closed-form antiderivatives `int_odd_model`, `int_even_model` verified against
  `mpmath` to $\le$ 9e-16 relative for $y$ up to $10^5$.
- `cumulative_theta` converges at 2nd order (1.26e-5 → 3.06e-9 over 64× refinement).
  Adequate for a first $\lambda$; will need replacing for CAP-grade seeds.

### Known remaining obstacle

The cumulative integral for $U$ has an integrand $\sim(\pi-\theta)^{p-2}$, which is
**not** integrable at the endpoint. One subtracted far-field term is not enough; either
a multi-term asymptotic expansion or the $U=HF$ route (one Hilbert transform, no
cumulative quadrature) is required. The $U=HF$ route needs $H$ of a *growing* function,
so the Cayley identity as written does not apply directly. This is the next design
decision.

## Files

| File | Role |
|---|---|
| `hilbert.py` | Cayley-map Hilbert transform, grid construction. |
| `ccf_profile.py` | Profile equation, far-field models, residual, solver. |
| `tests/test_hilbert.py` | Exact-pair and $H^2=-I$ validation. |
| `tests/test_components.py` | Quadrature, closed forms, $H\Omega$ vs `mpmath` PV integrals, self-convergence. |
