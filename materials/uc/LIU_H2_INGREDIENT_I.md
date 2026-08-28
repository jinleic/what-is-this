# Liu H2 ingredient (i): audit of the smooth A/C modes

## Status legend and executive result

**PROVED.** Every audit claim below begins with one of the required labels. A label applies to the whole paragraph, bullet, displayed formula, or table row that it introduces. Text inside source-code and verbatim-output blocks is evidence and retains the spelling used by the source.

**PROVED.** On the exact active-mean smooth chart, the apparent endpoint degeneracy is a coordinate degeneracy: both the quadratic objective excess and the quadratic tube distance contain the same factor \(q(1-q)\). The factor cancels exactly between the two quadratic forms for every \(0<q<1\); at \(q=0,1\), the corresponding \(d\)-direction is an inactive-component gauge.

**REFUTED.** The number `0.7004675091530868...` printed by `liu9_tube.py` is not the second-order ceiling for the denominator-cleared gap `numerator-EHX` requested in the tube lemma. It is the ceiling for the normalized excess \(\Phi-1\). The raw-gap ceiling has the additional point factor \(H_*=p h(x)\); the independently observed split-mode limit is `0.3871250877692325155...`.

**MACHINE VERIFIED.** A separate 320-bit Arb run encloses the exact conversion factor as

\[
 H_*=EHX(P_*)=p h(x)\in
 [0.552666730020487676432677600610273826374730816489206495137844711108831
 \mathbin{+/-}4.15\times10^{-70}].
\]

**MACHINE VERIFIED.** The same 320-bit calculation proves that the split ratio is the smaller of the two normalized ratios and encloses its raw-gap conversion as

\[
 k_G^{(2)}\in
 [0.3871250877692325155318291285580680847070641467880468753578147880128
 \mathbin{+/-}5.24\times10^{-68}].
\]

**PROVED.** Positivity of the point coefficients \(A\) and \(C\), together with an analytic factorization uniform in \(q\), proves a qualitative existential smooth-chart lemma: some neighborhood and some positive raw-gap constant exist on that chart.

**OPEN.** No current source supplies an explicit certified chart radius, a one-sided uniform remainder bound, or an explicit raw-gap \(\kappa\) on that neighborhood. Those, rather than another point Hessian calculation, are the missing machine-checkable content of ingredient (i).

**PROVED.** The primary classification is **analytic**. The core work is a uniform Taylor/divided-difference estimate for entropy expressions. Quotient gauges make the statement geometric, and Arb can certify explicit constants, but there is no combinatorial obstruction.

## 1. Definitions and the smooth quotient chart

**PROVED.** Let

\[
 h(t)=-t\log t-(1-t)\log(1-t),\qquad
 \pi(y,z)=yz\bigl(1+(1-y)(1-z)\bigr).
\]

**PROVED.** For a nine-coordinate point, let \(N=(1-\beta)EHXY+\beta EHPI\), \(H=EHX\), \(\Phi=N/H\), and \(G=N-H\). Thus `gap` in the requested tube lemma is \(G\), not \(\Phi-1\).

**PROVED.** At every point where \(H>0\), not merely in the limit at \(P_*\), the identity

\[
 G=N-H=H(\Phi-1)
\]

**PROVED.** This identity is exact. At \(H=0\), \(G=N-H\) remains the continuous gap expression while \(\Phi=N/H\) is undefined, so the quotient identity is not a statement there. In the smooth chart below, \(H>0\).

**PROVED.** Put \(m=px\), \(H_*=p h(x)\), and use the exact active-mean chart

\[
\begin{aligned}
 a(s)&=m/s, &x_0&=s-qd, &x_1&=s+(1-q)d,\\
 \Psi(q,s,d,r)&=(a(s),r,q,x_0,0,0,x_1,0,0),
\end{aligned}
\]

**PROVED.** Here \(0\le q\le1\), \(0\le r\le1-a(s)\), and the supports lie in \([0,1]\). The third mass is \(1-a(s)-r\).

**PROVED.** This chart is mean-active exactly, because

\[
 (1-q)a(s)x_0+q a(s)x_1=a(s)s=m.
\]

**PROVED.** The coordinate \(r\) only repartitions two atoms at support zero. At \(d=0\), changing \(q\) mixes two identical laws. At \(q=0\), changing \(d\) changes only inactive \(P_1\); at \(q=1\), it changes only inactive \(P_0\). These are the stated gauges.

**PROVED.** Write \(u=s-x\). The natural second-order quotient seminorm is

\[
 \|(u,d)\|_{q,*}^{2}
   =p x^2 u^2+(p^2+p x^2)q(1-q)d^2.
\]

**PROVED.** It is a genuine norm after quotienting the endpoint inactive-law directions and is exactly the quadratic part of the squared tube distance on this chart.

## 2. Reproduction ledger

**MACHINE VERIFIED.** All four required scripts were run exactly once from the repository root, sequentially, with one-core variables set to `1`. The tube script used `--skip-complement` because only its local-analysis block was under audit. Timing is `/usr/bin/time -p` wall time (`real`).

```sh
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export BLIS_NUM_THREADS=1

/usr/bin/time -p ./.venv/bin/python -I -B uc/liu9_tube.py --skip-complement
/usr/bin/time -p ./.venv/bin/python -I -B uc/liu9_binding.py
/usr/bin/time -p ./.venv/bin/python -I -B uc/liu9_reduce.py
/usr/bin/time -p ./.venv/bin/python -I -B uc/liu9_endpoint.py
```

| Status | Script | Exact arguments | Verbatim `real` time |
|---|---|---|---:|
| MACHINE VERIFIED | `uc/liu9_tube.py` | `--skip-complement` | `0.39` s |
| MACHINE VERIFIED | `uc/liu9_binding.py` | none | `2.05` s |
| MACHINE VERIFIED | `uc/liu9_reduce.py` | none | `2.36` s |
| MACHINE VERIFIED | `uc/liu9_endpoint.py` | none | `64.53` s |

**MACHINE VERIFIED.** The timing trailers were, verbatim:

```text
# liu9_tube.py --skip-complement
real 0.39
user 0.33
sys 0.05

# liu9_binding.py
real 2.05
user 1.99
sys 0.04

# liu9_reduce.py
real 2.36
user 2.30
sys 0.04

# liu9_endpoint.py
real 64.53
user 64.30
sys 0.14
```

**MACHINE VERIFIED.** Two additional unit probes were run once each under the same one-core environment:

```sh
/usr/bin/time -p ./.venv/bin/python -I -B /tmp/liu_hstar_arb.py
/usr/bin/time -p ./.venv/bin/python -I -B /tmp/liu_smooth_units_arb.py
```

**MACHINE VERIFIED.** Each probe took `real 0.37` seconds. The first set `ctx.prec = 320`, called `solve_equation_parameters(110)` and `certify_equation_parameters`, and evaluated `parameters.p * h_arb(parameters.x)`. The second also called `certify_curvatures`, reconstructed the two source ratios, proved `split_ratio < face_ratio`, and multiplied the smaller ratio by \(H_*\). Their relevant output was, verbatim:

```text
MACHINE VERIFIED [320-bit Arb]: H*=EHX(P*)=p*h(x) in [0.552666730020487676432677600610273826374730816489206495137844711108831 +/- 4.15e-70]
MACHINE VERIFIED [320-bit Arb]: normalized face ratio in [0.870199021581135103153079954399866437789814366374645938014244017910 +/- 6.52e-67]
MACHINE VERIFIED [320-bit Arb]: normalized split/common ratio in [0.7004675091530868229623611505681430164001846592010210319298715597306 +/- 5.02e-68]
MACHINE VERIFIED [320-bit Arb]: raw split/common ratio in [0.3871250877692325155318291285580680847070641467880468753578147880128 +/- 5.24e-68]
```

## 3. What the source chain actually proves

### 3.1 `liu9_tube.py`: point ratios, not a tube theorem

**PROVED.** `measure_local_attempts` computes two pointwise quadratic ratios. The source at lines 321--333 is:

```python
    # On the active-mean mode, dist^2=p*x^2*ds^2+o(ds^2).  On the paired
    # component mode, dist^2=q(1-q)*(p^2+p*x^2)*dd^2+o(dd^2).
    # Since the objective expansion uses one half of the Hessian, these are
    # the two exact quadratic comparison ratios for the displayed smooth modes.
    face_ratio = curvatures.face / (2 * parameters.p * parameters.x**2)
    split_ratio = curvatures.split_unit / (
        2 * (parameters.p**2 + parameters.p * parameters.x**2)
    )
    ratio_lower = min(face_ratio.lower(), split_ratio.lower(), key=float)
    ratio_upper = min(face_ratio.upper(), split_ratio.upper(), key=float)
    smooth_kappa = arb_hull(ratio_lower, ratio_upper)
    if not smooth_kappa > arb_fraction(Fraction(3, 5)):
        raise AssertionError("smooth quadratic comparison did not certify kappa>3/5")
```

**PROVED.** The comment correctly says “objective expansion.” Since `curvatures.face` and `curvatures.split_unit` are Hessian coefficients of \(\Phi\), these ratios compare \(\Phi-1\) with `dist^2`. They do not compare \(G=H(\Phi-1)\) with `dist^2`.

**PROVED.** `print_local_analysis` merely reports the imported point calculation. Its relevant source at lines 612--624 is:

```python
    curvatures = certify_curvatures(parameters)
    print("2. LOCAL TUBE ATTEMPT")
    print(f"PROVED [Arb, imported binding calculation]: A in {curvatures.face}")
    print(f"PROVED [Arb, imported binding calculation]: C in {curvatures.split_unit}")
    print("PROVED [exact reduced coordinates]: the Hessian is")
    print("   diag(0,A,C*q*(1-q),0) in (q,s,d,r); the q and r nulls are gauges.")
    print("PROVED [distance comparison on the two smooth displayed modes]:")
    print("   face ratio=A/(2*p*x^2), split ratio=C/[2*(p^2+p*x^2)], so")
    print(f"   their common quadratic kappa ceiling is in {attempts[0].smooth_kappa} > 3/5.")
    print("PROVED [liu9_endpoint.py, prior certificate]: D(Q)>=0 at q=0 and q=1,")
    print("   with equality exactly at Q=P* modulo its coordinate gauges.  This handles")
    print("   the sign of the endpoint first variation, but does not provide a uniform")
    print("   positive lower coefficient or a q-uniform higher-order remainder.")
```

**MACHINE VERIFIED.** The corresponding local output from the run was, verbatim:

```text
2. LOCAL TUBE ATTEMPT
PROVED [Arb, imported binding calculation]: A in [0.742135145133455480151420280758908219706652371733399528559468575887 +/- 4.65e-67]
PROVED [Arb, imported binding calculation]: C in [1.7160696505225689309702860986235592314301103639931179041197173441407 +/- 8.41e-68]
PROVED [exact reduced coordinates]: the Hessian is
   diag(0,A,C*q*(1-q),0) in (q,s,d,r); the q and r nulls are gauges.
PROVED [distance comparison on the two smooth displayed modes]:
   face ratio=A/(2*p*x^2), split ratio=C/[2*(p^2+p*x^2)], so
   their common quadratic kappa ceiling is in [0.7004675091530868229623611505681430164001846592010210319298715597306 +/- 5.02e-68] > 3/5.
PROVED [liu9_endpoint.py, prior certificate]: D(Q)>=0 at q=0 and q=1,
   with equality exactly at Q=P* modulo its coordinate gauges.  This handles
   the sign of the endpoint first variation, but does not provide a uniform
   positive lower coefficient or a q-uniform higher-order remainder.
```

**OPEN.** No line in this block chooses a neighborhood radius, bounds a remainder, or sets `tube_certified=True`; the same run explicitly printed `tube inequality certified=NO` for every requested radius.

**PROVED.** The units of the quantities involved are:

| Status | Quantity | Units relative to `dist^2` |
|---|---|---|
| PROVED | `face_ratio=A/(2*p*x^2)` | \((\Phi-1)/\mathrm{dist}^2\) at second order |
| PROVED | `split_ratio=C/[2*(p^2+p*x^2)]` | \((\Phi-1)/\mathrm{dist}^2\) at second order |
| PROVED | `smooth_kappa` and `smooth_kappa_ceiling` | the minimum of those normalized-objective units |
| PROVED | \(H_*\,\texttt{face_ratio}\), \(H_*\,\texttt{split_ratio}\) | raw \(G/\mathrm{dist}^2\) at second order |
| PROVED | `boundary_barrier` and a boundary-layer coefficient derived from a change in `gap` | raw denominator-cleared-gap units |
| PROVED | endpoint \(D(Q)\) | normalized objective first variation after active-mean compensation; raw-gap first variation at equality is \(H_*D(Q)\) |

**PROVED.** The factor converting point quadratic coefficients is exactly \(p h(x)\), not a rounded decimal. Away from the equality point the exact identity uses the point-dependent \(H(P)\), so multiplication by the fixed \(H_*\) is a second-order conversion at \(P_*\), not an exact finite-displacement replacement for \(H(P)\).

**PROVED.** The zero-boundary layer's admissible \(\kappa\), being derived from a change in denominator-cleared `gap`, is in raw-gap units. Comparing its approximately `0.873` ceiling with the old normalized `0.700...` smooth number mixes units; the corrected comparison is with raw `0.387125...`, so the conclusion that the boundary layer is non-binding survives with greater margin.

**COMPUTATIONAL EVIDENCE.** The lead agent reports the same unit correction for the mirror-layer \(\kappa\) column: against raw `0.387125...`, every certified `(rho,t0)` case is non-binding except the marginal `(0.1,1/64)` case. That file and its certificate remain under lead-agent ownership and are not re-audited here.

### 3.2 `liu9_binding.py`: exact contents of `CurvatureCertificate`

**PROVED.** The data container at lines 95--103 has only seven Arb-valued fields:

```python
@dataclass(frozen=True)
class CurvatureCertificate:
    face: arb
    split_unit: arb
    boundary_barrier: arb
    gain: arb
    loss: arb
    loss_over_gain: arb
    gain_minus_loss: arb
```

**PROVED.** `certify_curvatures`, the only definition of that name in `uc/`, is at lines 343--372:

```python
def certify_curvatures(parameters: ArbParameters) -> CurvatureCertificate:
    x, p, beta = parameters.x, parameters.p, parameters.beta
    face = parameters.mean * parameters.core_second

    x_jet = _SecondJet(x, 1, 0)
    diagonal_protocol = x_jet * x_jet * (
        1 + (1 - x_jet) * (1 - x_jet))
    protocol_diagonal_second = _h_jet(diagonal_protocol).second
    f_xx = _hpp_arb(x * x) * x**2
    split_unit = (
        2 * (1 - beta) * p * f_xx
        + beta * p * protocol_diagonal_second
        - _hpp_arb(x)
    ) / _h_arb(x)

    boundary_barrier = 2 * p * x * (1 + beta * (1 - x)) - 1
    gain = p * (2 * p * _f_xx_arb(x, x, beta) - _hpp_arb(x))
    loss = -2 * beta * p**2 * _k_xy_arb(x, x)
    if not (face > 0 and split_unit > 0 and boundary_barrier > 0
            and gain > 0 and gain - loss > 0):
        raise AssertionError("Arb did not certify the required positive curvatures")
    return CurvatureCertificate(
        face=face,
        split_unit=split_unit,
        boundary_barrier=boundary_barrier,
        gain=gain,
        loss=loss,
        loss_over_gain=loss / gain,
        gain_minus_loss=gain - loss,
    )
```

**PROVED.** If

\[
 \operatorname{core}(z)=
 \frac{(1-\beta)h(z^2)+\beta h(\pi(z,z))}{z h(z)},
\]

**PROVED.** The diagonal two-atom objective on the active-mean chart is \(m\operatorname{core}(s)\), so `face = mean * core_second` is \(A=\partial_s^2\Phi\) at \(s=x\).

**PROVED.** With \(F(y,z)=(1-\beta)h(yz)+\beta h(\pi(y,z))\), the code defines

\[
 \mathrm{GAIN}=p\{2pF_{11}(x,x)-h''(x)\},\qquad
 \mathrm{LOSS}=-2\beta p^2 k_{12}(x,x).
\]

**PROVED.** Symmetry of \(k\) and the diagonal second derivative give the exact identity

\[
 C=\texttt{split\_unit}
   =\frac{\mathrm{GAIN}-\mathrm{LOSS}}{H_*}.
\]

**PROVED.** `boundary_barrier` is the point coefficient

\[
 B=2px\bigl(1+\beta(1-x)\bigr)-1
\]

**PROVED.** This is the leading denominator-cleared \(y\log(1/y)\) boundary coefficient. The function asserts \(A>0\), \(C>0\), \(B>0\), `gain > 0`, and `gain-loss > 0`. It computes, but does not separately assert in its `if`, `loss > 0` or `loss/gain < 0.212`.

**MACHINE VERIFIED.** The binding run printed these relevant lines verbatim:

```text
PROVED [exact active-mean coordinates]: write
   x0=s-q*d, x1=s+(1-q)*d, p=m/s, and let r split the zero mass.
   Then the mean is identically m.  At (s,d)=(x,0), in coordinates
   (q,s,d,r), the reduced gradient is zero and the Hessian is
   diag(0,A,C*q*(1-q),0).
   The q null direction is exact P0=P1 invariance; the r null direction
   merely repartitions coincident zero atoms.
NUMERICAL [mpmath 90 dps derivative check]: |reduced gradient|=1.4111115e-100
PROVED [Arb eigenvalue coefficient]: A in [0.742135145133455480151420280758908219706652371733399528559468575887 +/- 4.65e-67]
PROVED [Arb eigenvalue coefficient]: C in [1.7160696505225689309702860986235592314301103639931179041197173441407 +/- 8.41e-68]
PROVED [Arb stationarity enclosure]: core'(x) in [+/- 7.79e-69]
```

```text
PROVED [one-sided positive-mass box directions]: every zero-support atom
   activated into the domain contributes, to the denominator-cleared gap,
   a leading positive multiple of y*log(1/y).  Its common coefficient is
   B=2*p*x*[1+beta*(1-x)]-1 in [0.27277669042050611463956575037617414476920817359978819894057041977347 +/- 4.08e-69].
```

```text
PROVED [Arb, the only moving interior atom, unit displacement]:
   GAIN in [1.2028972945610509121031416558621567602779959830547319954976714812191 +/- 6.09e-68]
   LOSS in [0.2544826923193416700883533887578149288295558545367558606436881845036 +/- 2.52e-68]
   LOSS/GAIN in [0.2115581217698265043287766506788360449391645167172127099168001457450 +/- 6.09e-68]
   GAIN-LOSS in [0.9484146022417092420147882671043418314484401285179761348539832967155 +/- 6.17e-68]
PROVED [comparison with the earlier non-binding probe]: LOSS/GAIN is <0.212
   here, not >4.2848.  The paired curvature is
   Phi''(0)=q*(1-q)*C >0 for every 0<q<1.
```

**MACHINE VERIFIED.** For the supplied 320-bit `ArbParameters`, every returned interval encloses the displayed formula and the asserted sign tests passed. This is exactly what `CurvatureCertificate` guarantees operationally.

**OPEN.** `CurvatureCertificate` has no radius, chart domain, norm, entropy lower bound, Hessian-variation bound, Taylor remainder, distance remainder, or coverage statement. It therefore does not by itself guarantee any inequality on a neighborhood.

**PROVED.** The printed Hessian identity is an algebraic statement encoded as output text; `print_reduced_hessians` numerically checks only the reduced gradient. The interval formulas for \(A\) and \(C\) are consistent with the identity, but the source does not run a symbolic Hessian comparison against the nine-variable objective.

### 3.3 `liu9_reduce.py`: corrected attribution and the genuine \(q\)-factor

**PROVED.** The literal `(q,s,d,r)` coordinate block and `diag(0,A,C*q*(1-q),0)` do not occur in `liu9_reduce.py`; they occur in `liu9_binding.py` and are reprinted by `liu9_tube.py`. `liu9_reduce.py` instead supplies two useful structural identities.

**PROVED.** First, `_split_terms` at lines 193--201 gives the exact mixture decomposition:

```python
    numerator0 = (one - beta) * j_mu + beta * k_mu
    phi0 = numerator0 / entropy
    correction = beta * qbar * q * k_nu / entropy
    objective = phi0 + correction
    gap = numerator0 + beta * qbar * q * k_nu - entropy

    # This is the exact residual of the literal prior claim
    # Phi=Phi0+beta*qbar*q*R(nu).
    claimed_r_residual = beta * qbar * q * (k_nu / entropy - r_nu)
```

**MACHINE VERIFIED.** The run printed the exact identity verbatim:

```text
PROVED [direct expansion of Liu (81)]: with H(mu)=int h(x)dmu,
   J(mu)=iint h(xy)dmu dmu, K(x,y)=h(xy[1+(1-x)(1-y)]),
   mu=(1-q)P0+qP1 and nu=P1-P0,
   Phi = Phi0(mu) + beta*q*(1-q)*K(nu)/H(mu),
   Phi0(mu)=[(1-beta)J(mu)+beta*K(mu)]/H(mu).
PROVED [same expansion, denominator-free]:
   G:=numerator-H = G0(mu)+beta*q*(1-q)*K(nu).
PROVED [degenerate fibre]: if P0=P1 then nu=0, so the correction
   vanishes and the source objective is invariant in q for every q in [0,1].
```

**PROVED.** Second, the script states the true in-domain paired-support path with the same weighted displacements used by the \(d\)-chart.

**MACHINE VERIFIED.** Its run printed, verbatim:

```text
PROVED [true in-domain paired-atom path]: for a three-atom P=sum ai*delta_bi,
   set P0(e)=sum ai*delta_(bi-q*e*di) and
       P1(e)=sum ai*delta_(bi+(1-q)*e*di).
   This remains inside the nine variables for all sufficiently small |e|.
   Although its coordinate mean is fixed, its mixture measure changes at e^2.
PROVED [explicit in-domain second variation]: put
   F(x,y)=(1-beta)h(xy)+beta*K(x,y), Phi0=F(P,P)/H(P), and
   G_atom=sum_i ai*di^2[2 sum_j aj*F_xx(bi,bj)-Phi0*h''(bi)].
   For <dot_nu,f>=sum_i ai*di*f'(bi),
   Phi''(0)=q*(1-q)/H(P) * [G_atom+2*beta*K(dot_nu)]
           =q*(1-q)/H(P) * [GAIN-LOSS],
   LOSS=2*beta*(r-D), r=-R(dot_nu).
```

**PROVED.** At \(P_*\), only the atom at \(x\) moves in the smooth split mode; the atoms at zero have displacement zero. The last formula therefore independently explains why the objective Hessian has the factor \(q(1-q)\).

### 3.4 `liu9_endpoint.py`: exact scope of \(D(Q)\ge0\)

**PROVED.** At \(q=0\), the active law is \(P_*\) and the inactive law is an admissible three-atom law \(Q\) sharing one of the two mass representations of \(P_*\); \(q=1\) is symmetric. After eliminating the active-law mean tangent, the endpoint first-variation functional is

\[
 D(Q)=\frac{2(1-\beta)[J(P_*,Q)-J(P_*,P_*)]
 +\beta[K(Q,Q)-K(P_*,P_*)]-[H(Q)-H_*]}{H_*}
 -\operatorname{core}[M(Q)-m].
\]

**MACHINE VERIFIED.** The endpoint run printed its definition and domain verbatim:

```text
PROVED [definitions and active-mean tangent elimination]:
   h(z)=-z*log(z)-(1-z)*log(1-z), with h(0)=h(1)=0;
   H(Q)=integral h(y)dQ(y), M(Q)=integral y dQ(y),
   J(P,Q)=double_integral h(y*z)dP(y)dQ(z), and
   K(P,Q)=double_integral h(y*z*[1+(1-y)(1-z)])dP(y)dQ(z).
   P*=p*delta_x+(1-p)*delta_0, m=p*x, H*=p*h(x), where
   x^4-2*x^3+3*x^2-1=0 and p=h(x)/h(x^2).
   Since x^2[1+(1-x)^2]=1-x^2, core=1/(p*x) exactly.
   For every admissible inactive law Q,
   D(Q)={2(1-beta)[J(P*,Q)-J(P*,P*)]
         +beta[K(Q,Q)-K(P*,P*)]-[H(Q)-H*]}/H*
         -core*[M(Q)-m].
PROVED [endpoint domain, up to atom permutations]:
   Q=p*delta_u+r*delta_v+(1-p-r)*delta_w, 0<=r<=1-p,
   or Q=r*delta_u+(p-r)*delta_v+(1-p)*delta_w, 0<=r<=p,
   with (u,v,w) in [0,1]^3.  Each component is four-dimensional.
   At q=0 this is the arbitrary inactive law P1 while P0=P*; at q=1
   it is the arbitrary inactive law P0 while P1=P*.  The two cases are
   identical by exchanging P0,P1.  The -core term is exactly the active
   mean compensation; no unconstrained mean direction was discarded.
```

**MACHINE VERIFIED.** The completed default certificate printed, verbatim:

```text
PROVED [local expansion plus rigorous Arb B&B]: D(Q)>=0 on both
   four-dimensional endpoint domains, at q=0 and q=1.  Equality holds
   exactly when the inactive law Q equals P* (including its coordinate
   gauges and atom permutations).
PROVED [certificate totals]: evaluated=406054; local-cleared=609; infeasible=74699; residual=0; weakest complement margin=[1.9397933871402502843008343056827144710098812137117402781589039021479817738851792363169400436950e-8 +/- 7.00e-104]; wall=64.148 seconds.
```

**PROVED.** This theorem is an endpoint sign and equality theorem for the first variation. It does not state a constant \(c_D>0\) such that \(D(Q)\ge c_D\delta(Q)^2\), so it cannot alone be substituted into `gap >= kappa * dist^2`.

**CONDITIONAL.** Once local coercivity at every equality representation is quantitative, strict complement margins and compactness reduce such a \(D/\delta^2\) bound to finite bookkeeping. Thus a separate loss of a factor \(q\) is not the obstruction; extraction of an explicit ratio for inactive laws outside the smooth chart remains a certification task.

## 4. Independent structural derivation

### 4.1 Objective factorization

**PROVED.** Let \(F(q,s,d)=\Phi(\Psi(q,s,d))\) and \(f(s)=F(q,s,0)\). The function \(f\) is independent of \(q\), and

\[
 f(s)=m\,\operatorname{core}(s),\qquad f(x)=1,
 \qquad f'(x)=0,
 \qquad f''(x)=A.
\]

**PROVED.** On a fixed support box around \(x\), all nonzero entropy arguments stay in the interior of \((0,1)\), so \(F\) is real analytic in \((q,s,d)\). Component exchange gives \(\partial_dF(q,s,0)=0\) for every \(q,s\). Moreover, `q=0` and `q=1` make the displaced component inactive, so \(F(0,s,d)=F(1,s,d)=f(s)\). Analytic divisibility therefore gives

\[
 F(q,s,d)=f(s)+q(1-q)d^2 R(q,s,d),
\]

**PROVED.** The function \(R\) extends continuously, indeed analytically, to \(q\in\{0,1\}\) and \(d=0\).

**PROVED.** The point Hessian identity implies \(R(q,x,0)=C/2\) for every \(q\). Consequently, uniformly for \(q\in[0,1]\),

\[
 \Phi(\Psi)-1
 =\frac{A}{2}u^2+\frac{C}{2}q(1-q)d^2
 +O\!\left((|u|+|d|)\,[u^2+q(1-q)d^2]\right).
\]

**PROVED.** The weighted form of the remainder, rather than an unweighted `O(d^3)`, is the key endpoint-uniform statement.

### 4.2 Exact distance and its quadratic part

**PROVED.** On the same chart, with \(a=m/s\), the distance is exactly

\[
\begin{aligned}
 \operatorname{dist}^2={}&(1-q)\left[a^2q^2d^2
       +a\{x_0(x_0-x)\}^2\right]\\
 &+q\left[a^2(1-q)^2d^2
       +a\{x_1(x_1-x)\}^2\right].
\end{aligned}
\]

**PROVED.** Expanding at \((s,d)=(x,0)\) gives, uniformly in \(q\),

\[
 \operatorname{dist}^2
 =p x^2u^2+(p^2+p x^2)q(1-q)d^2
 +O\!\left((|u|+|d|)\,[u^2+q(1-q)d^2]\right).
\]

**PROVED.** Thus the objective and distance have identical \(q(1-q)\) weights in the split mode.

**PROVED.** For the pure \(d\) path \(u=0\), both the full finite-displacement objective excess and the full finite-displacement distance are exactly divisible by \(q(1-q)d^2\). Their residual quotients can still depend on \(q\) and \(d\). What is exactly \(q\)-independent is the ratio of the two quadratic coefficients at \(d=0\); the numerical finite-\(d\) ratios approach that value asymptotically as \(d\to0\).

### 4.3 Normalized objective versus raw gap

**PROVED.** Since \(G=H(\Phi-1)\), \(H(\Psi)=H_*+O(|u|+|d|)\), and \(H_*>0\),

\[
 G(\Psi)
 =\frac{H_*A}{2}u^2
  +\frac{H_*C}{2}q(1-q)d^2
  +O\!\left((|u|+|d|)\,[u^2+q(1-q)d^2]\right).
\]

**PROVED.** The two point ratios are therefore

\[
 k_{\Phi}^{(2)}=
 \min\left\{\frac{A}{2px^2},
             \frac{C}{2(p^2+px^2)}\right\},
 \qquad
 k_G^{(2)}=p h(x)\, k_{\Phi}^{(2)}.
\]

**PROVED.** The first ratio is what `liu9_tube.py` prints; the second is the ceiling relevant to the requested raw-gap lemma. The conversion factor in this equality is the exact expression \(p h(x)\), enclosed by Arb above.

## 5. Independent 160-decimal-digit experiment

**COMPUTATIONAL EVIDENCE.** The experiment used the exact equation-defined parameters, the repository's independent objective transcription, and the repository's distance function. For each \(q\), it fixed \(s=x\), used \(x_0=x-qd\), \(x_1=x+(1-q)d\), and set the two zero-mass representatives to support zero. Every row is mean-feasible up to the displayed 170-digit solve residual.

**COMPUTATIONAL EVIDENCE.** The exact execution command was:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
/usr/bin/time -p ./.venv/bin/python -I -B /tmp/liu_h2_ingredient_i_experiment.py
```

**COMPUTATIONAL EVIDENCE.** The temporary script executed by that command was exactly:

```python
import os
import sys

for name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ[name] = "1"

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/math/uc")

import mpmath as mp

from liu9_binding import solve_equation_parameters
from liu9_objective import evaluate_mpmath, mean_value
from liu9_tube import distance_squared_mp

mp.mp.dps = 180
parameters = solve_equation_parameters(170)
x = parameters.x
p = parameters.p
beta = parameters.beta
r = (1 - p) / 2


def h(z):
    if z == 0 or z == 1:
        return mp.mpf(0)
    return -z * mp.log(z) - (1 - z) * mp.log(1 - z)


def hpp(z):
    return -1 / (z * (1 - z))


def pi_diag(z):
    return z * z * (1 + (1 - z) ** 2)


protocol_diagonal_second = mp.diff(lambda z: h(pi_diag(z)), x, 2)
f_xx = hpp(x * x) * x * x
C = (
    2 * (1 - beta) * p * f_xx
    + beta * p * protocol_diagonal_second
    - hpp(x)
) / h(x)
Hstar = p * h(x)
distance_coefficient = p * p + p * x * x
objective_limit = C / (2 * distance_coefficient)
gap_limit = Hstar * objective_limit

print("COMPUTATIONAL EVIDENCE [mpmath 160 dps pure-d experiment]")
print("C=" + mp.nstr(C, 40))
print("Hstar=" + mp.nstr(Hstar, 40))
print("objective_limit=" + mp.nstr(objective_limit, 40))
print("gap_limit=" + mp.nstr(gap_limit, 40))
print("q | d | mean-m | gap/dist^2 | (Phi-1)/dist^2")
for q_text in ("1e-1", "1e-3", "1e-6"):
    q = mp.mpf(q_text)
    for d_text in ("1e-2", "1e-4", "1e-6", "1e-8"):
        d = mp.mpf(d_text)
        x0 = x - q * d
        x1 = x + (1 - q) * d
        values = (p, r, q, x0, mp.mpf(0), mp.mpf(0), x1, mp.mpf(0), mp.mpf(0))
        terms = evaluate_mpmath(values, beta, dps=160)
        gap = terms.numerator - terms.ehx
        dist2 = distance_squared_mp(values, parameters)
        mean_residual = mean_value(values) - parameters.mean
        print(
            q_text
            + " | "
            + d_text
            + " | "
            + mp.nstr(mean_residual, 5)
            + " | "
            + mp.nstr(gap / dist2, 30)
            + " | "
            + mp.nstr((terms.objective - 1) / dist2, 30)
        )
```

**COMPUTATIONAL EVIDENCE.** The run first printed these predicted limits and took `real 0.39` seconds:

```text
C=1.71606965052256893097028609862355923143
Hstar=0.5526667300204876764326776006102738263747
objective_limit=0.7004675091530868229623611505681430164002
gap_limit=0.3871250877692325155318291285580680847071
```

**COMPUTATIONAL EVIDENCE.** The complete numerical table was:

| \(q\) | \(d\) | `mean-m` | `gap/dist^2` | `(Phi-1)/dist^2` |
|---:|---:|---:|---:|---:|
| `1e-1` | `1e-2` | `6.3345e-173` | `0.390566463301781037012228433073` | `0.706718554589141273478442738866` |
| `1e-1` | `1e-4` | `6.3345e-173` | `0.38715917869749465796924443029` | `0.700529195962305539701050656687` |
| `1e-1` | `1e-6` | `6.3345e-173` | `0.3871254286472998396694015266` | `0.7004681259410725611294713063` |
| `1e-1` | `1e-8` | `6.3345e-173` | `0.387125091178010068289925939061` | `0.700467515320958671913144883299` |
| `1e-3` | `1e-2` | `6.3345e-173` | `0.391442332159060020488210574815` | `0.708279438227503890248256336241` |
| `1e-3` | `1e-4` | `6.3345e-173` | `0.387167618554620718870246309449` | `0.700544464753483868112494076699` |
| `1e-3` | `1e-6` | `6.3345e-173` | `0.387125513014856681697396116261` | `0.700468278596229554126238223338` |
| `1e-3` | `1e-8` | `6.3345e-173` | `0.387125092021682536115426496084` | `0.700467516847506967580677050255` |
| `1e-6` | `1e-2` | `6.3345e-173` | `0.391451207560961291245256048646` | `0.708295228293389222985858295253` |
| `1e-6` | `1e-4` | `6.3345e-173` | `0.387167703723968988899097146697` | `0.700544618833182001134208650994` |
| `1e-6` | `1e-6` | `6.3345e-173` | `0.387125513866202379554379488217` | `0.700468280136659236833024753629` |
| `1e-6` | `1e-8` | `6.3345e-173` | `0.387125092030195958325910559489` | `0.700467516862911227692263319128` |

**COMPUTATIONAL EVIDENCE.** As \(d\to0\), all three \(q\)-rows converge to the same normalized limit `0.7004675091530868...` and the same raw-gap limit `0.3871250877692325...`, including \(q=10^{-6}\). No extra factor of \(q\) appears in either ratio.

## 6. Precise proposition that ingredient (i) must mean

**OPEN.** The following explicit, uniform proposition is the missing certified version of ingredient (i), rather than the point-Hessian statement currently printed.

> **Smooth A/C proposition.** There exist explicitly certified positive rational constants \(\varepsilon,\kappa_{\rm sm},M_G,M_D\) such that, for every \(q\in[0,1]\), every \(s,d\) satisfying
> \[
> |s-x|\le\varepsilon,\qquad |d|\le\varepsilon,
> \qquad x_0=s-qd,\ x_1=s+(1-q)d\in[19/32,25/32],
> \]
> and every \(r\in[0,1-m/s]\), the point \(\Psi(q,s,d,r)\) obeys
> \[
> G(\Psi)\ge \kappa_{\rm sm}\operatorname{dist}(\Psi)^2.
> \]
> The same constants work for all \(q\), including \(q=0,1\), after quotienting the inactive-law and repeated-zero-atom gauges. In addition, with \(u=s-x\), the proof certifies the one-sided weighted remainders
> \[
> \begin{aligned}
> G(\Psi)&=\frac{H_*}{2}\{Au^2+Cq(1-q)d^2\}+R_G,\\
> \operatorname{dist}^2&=px^2u^2+(p^2+px^2)q(1-q)d^2+R_D,\\
> R_G&\ge-M_G(|u|+|d|)\{u^2+q(1-q)d^2\},\\
> |R_D|&\le M_D(|u|+|d|)\{u^2+q(1-q)d^2\}.
> \end{aligned}
> \]

**PROVED.** The unweighted coordinate box `max(|s-x|,|d|) <= epsilon` is essential. A neighborhood defined only by small `dist` permits an inactive law to remain far from \(P_*\) when \(q\) is tiny. Such points are not in the smooth chart and belong to the endpoint-\(D\) part of the piecewise proof.

**PROVED.** The proposition's norm is the weighted quotient norm \(\|\cdot\|_{q,*}\), while its neighborhood topology is the ordinary sup norm in \((s,d)\). Conflating these two roles is precisely what makes the endpoint factor look degenerate.

**PROVED.** If the displayed remainder bounds hold and \(\varepsilon\) is small enough that the lower quadratic coefficient remains positive and the distance upper coefficient is finite, then an explicit \(\kappa_{\rm sm}>0\) follows immediately. This is what the ingredient must deliver to the integrated `gap >= kappa * dist^2` lemma.

**CONDITIONAL.** In the full piecewise assembly, the global local radius is the minimum radius for which every point is assigned to this smooth chart, the proved zero boundary layer, the mirror layer, or an endpoint inactive-law chart. The final \(\kappa\) is the minimum of the corresponding raw-gap constants.

## 7. Point statement versus neighborhood statement

| Status | Statement |
|---|---|
| MACHINE VERIFIED | At \((s,d)=(x,0)\), \(A\) and \(C\) lie in the printed positive Arb intervals. |
| PROVED | At that point, the normalized Hessian is `diag(0,A,C*q*(1-q),0)` in the stated chart. |
| PROVED | The `q` and `r` point nulls are representation gauges; endpoint `d` is also an inactive-law gauge. |
| PROVED | The normalized second-order comparison ceiling is `0.7004675091530868...`; the raw-gap split ceiling is `0.3871250877692325...`. |
| PROVED | Analytic factorization and compactness in \(q\) imply existence of some uniform smooth-chart neighborhood and some positive constant. |
| OPEN | No explicit certified \(\varepsilon\) or one-sided weighted remainder constant is present in the current scripts. |
| OPEN | The current scripts do not prove that every point of a `dist <= rho` tube belongs to this smooth chart; q-small far inactive laws require endpoint treatment. |
| OPEN | The current scripts do not state the final raw-gap normalization when combining this mode with the other pieces. |

## 8. Weakest sufficient replacement and remaining proof work

**PROVED.** A full interval enclosure of every third derivative, or a full neighborhood Hessian spectrum, is stronger than necessary. The weakest sufficient replacement is a direct one-sided certificate for one positive constant:

\[
 G(\Psi(q,s,d,r))-\kappa_{\rm sm}\operatorname{dist}^2(\Psi(q,s,d,r))\ge0
\]

**PROVED.** This certificate is required on one explicit smooth chart box, uniformly in \(q\), with equality only on its quotient gauges. Any positive rational \(\kappa_{\rm sm}\) is sufficient; optimizing toward the point ceiling is unnecessary.

**PROVED.** An almost equally weak and more convenient certificate is the one-sided weighted remainder bound in the proposition. It uses the exact factorization by \(q(1-q)d^2\), so it never divides an ordinary interval by an interval containing zero.

**OPEN.** To machine-prove that replacement, the remaining work is:

1. **OPEN.** Encode the following analytic extension with integral divided differences, so its Arb enclosure remains finite at \(q=0,1\) and \(d=0\):
   \[
   [\Phi(q,s,d)-\Phi(q,s,0)]/[q(1-q)d^2]
   \]
2. **OPEN.** Certify a positive lower bound for that extension and for the diagonal quotient `[Phi(q,s,0)-1]/(s-x)^2` on one explicit rational chart box.
3. **OPEN.** Certify a positive lower bound for \(H\) and an upper comparison between the exact distance and \(\|\cdot\|_{q,*}^2\); multiply by \(H\) before reporting a raw-gap \(\kappa\).
4. **OPEN.** Check the finite collection of zero-atom representations and atom permutations without introducing a second smooth convention, and handle mean slack by its positive constrained first variation.
5. **OPEN.** For q-small points whose inactive law is outside this chart, extract a quantitative \(D(Q)/\delta(Q)^2\) constant from the endpoint local certificates plus the already strict complement boxes. This is endpoint bookkeeping, not an unmatched \(q\)-power in the A/C mode.

**PROVED.** Therefore the answer to the structural question is: the \(q(1-q)\) factors genuinely cancel in the second-order ratio. The cancellation of the powers is exact at quadratic order for every \(q\), not merely an empirical asymptotic as \(q\to0\); only the finite-\(d\) residual quotient varies before its \(d\to0\) limit. Ingredient (iii) is not needed to repair a smooth-mode q-weight mismatch, although a quantitative endpoint ratio is still needed for inactive laws that are not close in the ordinary smooth-chart topology.
