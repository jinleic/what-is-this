# Pre-statement — Gate-B direct Lemma-D.4 restart with analytic unit-cubic anchors

Campaign: `20260831T082425Z_kg_direct_d4_restart`

This is the first file in this campaign directory.  It is committed before any
campaign Python invocation or numerical computation.  The stopped predecessor
`20260831T080542Z_kg_band_close2` is immutable and will not be modified.

## 1. Fixed target and conventions

The target is paper inequality (56)

`J(c,p)=E[(|p(Z)|-c|Z|)_+] <= d(c)`

for every unit cubic `p=sum_{j=0}^3 a_j psi_j`, with
`d(c)=(3 nu/2)(sqrt(1+c^2)-c)` and `nu=sqrt(2/pi)`.  The basis is the
orthonormal probabilist convention

`psi_0=1`, `psi_1=s`, `psi_2=(s^2-1)/sqrt(2)`,
`psi_3=(s^3-3s)/sqrt(6)`.

The certificate route is paper D.3: fold onto `[0,infinity)`, split `p=e+o`,
and use Lemma D.4 on

`f(x,y;t)=(x+y-t)_+ + (|x-y|-t)_+`, `x,y,t>=0`.

## 2. Independent analytic anchor, derived before code

The historical center is itself a unit cubic:

`p_*(s)=(4/5)psi_0(s)-(3/5)psi_2(s)=A-Bs^2`,

where

`A=4/5+3/(5sqrt(2))`, `B=3/(5sqrt(2))>0`.

For `c>0`, on the positive half-line the active set of
`(|p_*(s)|-cs)_+` is exactly `[0,r1] union [r2,infinity)`, with roots

`r1=(sqrt(c^2+4AB)-c)/(2B)`,
`r2=(sqrt(c^2+4AB)+c)/(2B)`.

This follows by solving `A-Bs^2-cs=0` before the sign change of `p_*`, and
`Bs^2-A-cs=0` after it.  With `P(r)=Phi(r)-1/2`, `Q(r)=1-Phi(r)`, and
`varphi=phi`, the exact closed form is

`J(c,p_*) = 2 [ (A-B)P(r1) + B r1 varphi(r1)
                  - c(varphi(0)-varphi(r1))
                + (B-A)Q(r2) + B r2 varphi(r2)
                  - c varphi(r2) ]`.

It uses only Arb `sqrt`, `exp`, and `erf` with outward endpoints.  A separate
80-digit mpmath quadrature split at the two roots is permitted only as
`COMPUTATIONAL-EVIDENCE`; it cannot establish any sign.

At 300 Arb bits, startup assertions are pre-anchored to the independently
supplied narrow windows:

- at `c=1.30`,
  `r1 in [0.7554755377609217537989341027304,
          0.7554755377609217537989341027306]`;
- `r2 in [3.8196049229026276928692596718514,
          3.8196049229026276928692596718516]`;
- `J(1.30,p_*) in
  [0.37483055458876357002453330452976048,
   0.37483055458876357002453330452976050]`;
- `d(1.45)-J(1.30,p_*) in
  [-0.002148854794213660928402532894742,
   -0.002148854794213660928402532894740]`;
- `d(1.30)-J(1.30,p_*) in
  [0.032236520487306532595117305408006,
   0.032236520487306532595117305408008]`.

The coarse c-pair `(1.30,1.45)` must therefore be rejected before any envelope
run: because the center belongs to the coefficient cell, no sound upper
envelope can have a positive `d(1.45)-U(1.30,B)` margin.  This is an analytic
control of semantics, not merely provenance.

## 3. Direct certified panel envelope

For a closed even-coefficient rectangle `B` and panel `P=[a,b] subset [0,8]`,
the instrument computes `E_P>=|e(s)|` by outward interval evaluation of the
affine even polynomial over `B x P`.  For every compatible unit cubic,

`|o(s)| <= W_P = sqrt(1-dist(0,B)^2) sqrt(K_o(b))`,

where `K_o=(5/2)s^2-s^4+s^6/6` and
`K_o'(s)=s((s^2-2)^2+1)>0` for `s>0`.

Lemma D.4 and `c_L s >= c_L a` give, directly and without an `mq` case split,

`f(|e|,|o|;c_L s)
 <= (E_P+W_P-c_L a)_+ + (|E_P-W_P|-c_L a)_+ = G_P`.

The certified upper bound is

`U(B,c_L)=sum_P [Phi(b)-Phi(a)] G_P + Tail_8`.

The exact CDF mass is already the inherited sound design.  The Gaussian
density decreases on `[0,infinity)`, so `phi(b)(b-a)` is a planted permissive
underweight, not an allowed alternative.

For every c-tile `[c_L,c_U]`, `J(c,p)<=J(c_L,p)` and
`d(c)>=d(c_U)`.  A leaf passes only under the strict separated-endpoint test

`d(c_U).lower() - U(B,c_L).upper() > 0`.

A band passes only if every even-disk leaf in every adjacent c-tile passes.
Failed internal parents do not enter the reported certified-leaf margin.

## 4. Mandatory startup controls and stop rule

The controls run in this order and abort on first failure:

1. analytic coefficient norm `((4/5)^2+(3/5)^2)=1` and `p_*` closed-form
   roots/J/margins in the windows above;
2. a root-partition check that the signs of the two quadratics at `0,r1,r2`
   establish exactly the stated active components;
3. Rule-17b root-cover count: exactly `132` boxes in the exact-rational
   `12 x 12` partition intersect the closed unit disk;
4. Gaussian-weight plant `phi(1)<Phi(1)-Phi(0)<phi(0)` and strict-margin
   positive/negative/overlap plants;
5. direct-envelope containment: on the radius-`0.0001` cell centered at
   `(0.8,-0.6)`, `U(B,1.30).upper() >= J(1.30,p_*).lower()`; a deliberately
   planted `J-10^-6` pseudo-upper-bound must be rejected;
6. direct-envelope equal-c anchor on that cell, c-pair `(1.30,1.30)`, `512`
   panels: strict positive margin required;
7. direct-envelope ratio-`1.02` tile anchor on that cell, c-pair
   `(1.30,1.326)`, `512` panels: strict positive margin required.

If either direct-envelope anchor fails, the campaign stops honestly and no
full band is run.  The impossible `(1.30,1.45)` pair is never used as a
positive anchor.

## 5. Fixed full sweep, ladders, and budget

Only after all startup controls pass, run these six full bands in the fixed
order:

1. `[4.083,6.0]`;
2. `[1.0,1.3]`;
3. `[1.30,1.45]`;
4. `[1.45,1.75]`;
5. `[1.75,3.5]`;
6. `[3.50,4.083]`.

Each band is tiled deterministically by adjacent ratio-`1.02` pairs
`[a,min(1.02a,c_hi)]`, with the final endpoint clamped exactly to `c_hi`.
There is no c-monotonicity shortcut across a whole named band.

The even disk uses the exact `12 x 12` root partition and binary longer-axis
splits, maximum depth `40`, minimum side `1/1024`.  The unique panel ladder is
the inherited start times powers of four, capped once at `32768`: start `2048`
for `c_L<=1`, `1024` for `1<c_L<=4`, `4096` for `4<c_L<=8`, `8192` above
`8`.  Precision is `256` Arb bits for the envelope and `300` bits for analytic
controls.  Each named band has a budget of `260000` distinct envelope
evaluations.  Budget exhaustion or one unresolved floor leaf earns only
`FAILURE TO CERTIFY`, with the exact frontier.  No post-result domain change is
allowed.

Every band artifact records total/pre-anchored root counts, tiles
attempted/certified/open, certified leaves, unresolved frontier, distinct
envelope evaluations, unique panel-count histogram, minimum certified-leaf
margin with Arb radius, and unresolved margin.  A startup assertion checks
that tile endpoints are adjacent and cover the exact named band without gaps.

## 6. Evidence and refutation gate

`MACHINE-VERIFIED` means outward Arb endpoints establish the stated sign over
the stated coefficient/c/panel/tail domain.  Float/mpmath values are
`COMPUTATIONAL-EVIDENCE`.  Paper Table-1 values are `REPORTED`.  Counts earn
machine-verified status only after the independent known-count startup assert.

An envelope margin `<=0` is only failure to certify.  Escalation to Main as a
candidate paper refutation requires an explicit unit cubic and c with a
separate certified lower bound `J.lower()>d.upper()`.  The analytic coarse
pair compares `J(1.30,p_*)` to `d(1.45)`, not to `d(1.30)`, so its negative
paired margin is not a violation of the pointwise paper inequality.

## 7. Pre-registered Rule-7 scope sentence

This campaign sweeps exactly all unit cubics in the orthonormal-probabilist
Hermite span `psi_0,...,psi_3` through the full even-coefficient disk, with all
compatible odd coefficients covered by `sqrt(1-dist(0,B)^2)`, over the six
full c-bands `[1.0,1.3]`, `[1.30,1.45]`, `[1.45,1.75]`, `[1.75,3.5]`,
`[3.50,4.083]`, and `[4.083,6.0]`, using adjacent ratio-`1.02` c-tiles, the
exact `12 x 12` root cover and stated binary radii down to side `1/1024`, the
stated panel ladders through `32768`, 256-bit Arb panels on `[0,8]`, and the
closed `[8,infinity)` tail; it separately evaluates only the explicit cubic
`p_*=(4/5)psi_0-(3/5)psi_2` at c `1.30` and d at c `1.30` and `1.45`, plus the
radius-`0.0001` anchor cell at c-pairs `(1.30,1.30)` and `(1.30,1.326)`, and it
does **not** search c in `[0.993405,1.0)`, c above `6.0`, other isolated cell
radii, non-unit cubics, degrees above three, the paper's C1/C3/splice
certificates, Gate-C septic/Krivine scheme families, or constructions outside
the unit-cubic threshold family.
