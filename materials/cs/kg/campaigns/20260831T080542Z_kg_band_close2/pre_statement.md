# Pre-statement — Gate-B certified band close-out, corrected paper-D.3 envelope

Campaign: `20260831T080542Z_kg_band_close2`

This file is the first file in the campaign directory and is committed before any
campaign Python invocation or numerical result.

## 1. Fixed target and primary-source route

The target is paper inequality (56): for every orthonormal-probabilist Hermite
unit cubic

`p(s) = sum_{j=0}^3 a_j psi_j(s)`, `sum a_j^2 = 1`,

certify

`J(c,p) = E[(|p(Z)| - c|Z|)_+] <= d(c)`,

where `d(c) = (3 nu/2)(sqrt(1+c^2)-c)` and `nu=sqrt(2/pi)`.
The conventions are `psi_0=1`, `psi_1=s`, `psi_2=(s^2-1)/sqrt(2)`, and
`psi_3=(s^3-3s)/sqrt(6)` (paper line 942).  The route is the paper's D.3
half-line fold (57), Lemma D.4 monotone envelope, Christoffel odd envelope,
interval panels on `[0,8]`, the stated closed tail majorant, and the band
pairing on paper line 1984.  These statements were re-read from the paper
before this campaign was created; no output shape was copied in their place.

For `s>=0`, write `p=e+o`, `x=|e(s)|`, `y=|o(s)|`, and `t=cs`.  The folded
integrand is exactly

`f(x,y;t)=(x+y-t)_+ + (|x-y|-t)_+`.

For an even-coefficient rectangle `B`, the instrument computes an outward Arb
upper bound `E_P >= |e(s)|` on every panel `P=[a,b]`.  Every compatible unit
cubic has

`|o(s)| <= W_P := sqrt(1-dist(0,B)^2) sqrt(K_o(b))`,

with `K_o(s)=(5/2)s^2-s^4+s^6/6`; its derivative is
`s((s^2-2)^2+1)>0`, so the right endpoint is the panel supremum.  Lemma D.4
and `t=c_L s >= c_L a` give the certified panel supremum

`G_P = (E_P+W_P-c_L a)_+ + (|E_P-W_P|-c_L a)_+`.

Thus

`J(c_L,p) <= sum_P (Phi(b)-Phi(a)) G_P + Tail_8`.

For a c-tile `[c_L,c_U]`, `J(c,p)<=J(c_L,p)` and
`d(c)>=d(c_U)` because both functions decrease in c.  A tile passes only if
the Arb upper bound on the left is strictly below the Arb lower bound on
`d(c_U)` for every even-disk leaf.  A full band passes only if every c-tile
and every disk leaf passes.  A single cell, including the startup anchor, can
never establish a band claim.

## 2. Source-audit corrections fixed before computation

The inherited frozen source and artifacts were read line by line.  This
campaign accepts none of their band verdicts by inheritance.

1. The inherited `try_cell` used `(d_lower-U).upper()>0`.  This is the
   permissive endpoint and is not a certificate of a positive margin.  This
   campaign uses `d.lower()-U.upper()>0` only.
2. The inherited ordered-`q` discussion treated
   `A2-r2` as a lower bound on `||e|^2-|o|^2|`; it is not: `A2` and `r2` are
   independent upper envelopes.  More importantly, no `mq=t^2` test is
   needed to use the tight expression above.  Lemma D.4 directly proves
   `f(x,y;t)<=f(E_P,W_P;t)`.  The inherited positive-branch code already
   equals this paper envelope; its fallback merely drops `-t` from the
   second hinge and is looser.  The pre-registered first change is therefore
   to delete the semantically unsupported `mq` branch and use the direct
   Lemma-D.4 envelope on every panel.  Panel subdivision at `mq=t^2` would
   refine a condition that is neither necessary nor established, so it will
   not be implemented unless this derivation is falsified.
3. Exact panel mass `Phi(b)-Phi(a)` is already present in the inherited code;
   it is not a new optimization in this campaign.  Also, on `[a,b]` with
   `0<=a<b`, Gaussian density decreases, so `phi(b)(b-a)` is a lower
   rectangle, not a supremum weight.  The sound rectangular comparison is
   `phi(a)(b-a)`.  Startup controls will plant `[0,1]` and require
   `phi(1) < Phi(1)-Phi(0) < phi(0)`; this rejects the permissive
   right-endpoint-weight variant.
4. The inherited successful-band margin accumulator included failed internal
   branch-and-bound parents.  This campaign reports minimum margins over
   certified leaves only and reports unresolved frontier margins separately.
5. The inherited frozen `b5_result.json` reports zero certified c-tiles and
   unresolved margins from about `-0.0283` to `-0.0194`, whereas the appended
   README text describes seven of eight tiles and `-2e-5` to `-5e-5`.  This
   campaign treats the frozen JSON as the only as-run evidence and will record
   the stale README wording inline rather than silently preserve it.

## 3. Startup anchors and counterfactual controls

All controls run before any band sweep and abort the campaign on failure.

* Basis anchor: the normalized reproducing kernel at zero has
  `(a0,a2)=(sqrt(2/3),-1/sqrt(3))`, norm exactly one and `p(0)^2=3/2`.
* Coverage-count anchor (Rule 17b): a `12 x 12` partition of `[-1,1]^2`
  contains exactly **132** boxes intersecting the closed unit disk.  This is
  fixed independently by integer distances
  `(5,4,3,2,1,0,0,1,2,3,4,5)/6` and is asserted at startup.
* Weight counterfactual: `[0,1]` must reject `phi(1)` as an upper rectangle and
  accept `phi(0)` as one, as specified in section 2.
* Margin counterfactuals: a planted upper bound strictly above `d`, and a ball
  straddling `d`, must both be rejected; a planted upper bound strictly below
  `d` must be accepted.
* Historical numerical anchor: before changing the copied source, the frozen
  inherited instrument is run on only the cell
  `(a0,a2)=(0.8,-0.6)+[-1e-4,1e-4]^2`, c-pair `(1.30,1.45)`, 512 panels.
  Its certified lower margin must lie in `[0.192,0.194]`, reproducing the
  previously recorded rounded `+0.1930`; otherwise the campaign stops.  This
  is a pipeline anchor only, not a band verdict.  The corrected instrument is
  then required to pass the same cell, with a margin no smaller than the old
  one because the replacement is pointwise no looser.

## 4. Fixed sweep, ladders, and budgets

The c-tiles use the pre-existing deterministic geometric ratio `1.02`, with
adjacent exact Arb endpoints and the final endpoint clamped to the named band.
The full even disk is covered by the 132-root partition above; unresolved
rectangles split along their longer axis to maximum depth 40 and minimum side
`1/1024`.  Every child retained is one whose closed rectangle intersects the
closed unit disk.

Unique panel ladders are the inherited starts multiplied by four and capped at
32768: start 2048 for `c_L<=1`, 1024 for `1<c_L<=4`, 4096 for
`4<c_L<=8`, and 8192 above 8.  Repeated evaluations at the 32768 cap are
removed.  Precision is 256 Arb bits.  The tail domain is exactly `[8,infinity)`
and uses the paper's closed majorant.  The per-band budget is 260,000 distinct
certified cell-envelope evaluations.  Budget exhaustion earns only
`FAILURE TO CERTIFY` plus the exact frontier.

Order, fixed before results:

1. startup controls and inherited single-cell anchor;
2. corrected single-cell anchor;
3. exact-CDF versus sound left-rectangle diagnostic on the anchor and the
   razor cell `(0.8,-0.6)+[-1e-4,1e-4]^2`, c-pair `(6.0,6.05)`, over the
   fixed ladder `n in {512,1024,2048,4096,8192}`;
4. full band `[4.083,6.0]`;
5. full bands `[1.0,1.3]`, `[1.30,1.45]`, `[1.45,1.75]`, `[1.75,3.5]`, and
   `[3.50,4.083]`, in that order.

Each full-band artifact must state c-tiles attempted/certified/open, certified
leaf count, unresolved frontier, distinct envelope evaluations, panel-count
histogram, worst certified leaf margin, and worst unresolved margin.  Counts
are semantic only after the startup coverage anchor passes.

## 5. Refutation trigger and evidence strength

An envelope upper bound exceeding `d` is only a failure to certify; it is not
evidence that the paper inequality is false.  Escalation to Main is triggered
only by an explicit unit cubic and c for which a separate certified *lower*
bound on `J(c,p)` is strictly greater than an Arb upper bound on `d(c)`, after
root-aware subdivision of every active interval.  No paper disagreement will
be recorded before that escalation.

`MACHINE-VERIFIED` is reserved for signs obtained from outward Arb endpoints
through the whole stated domain.  Float probes, timing, and heuristic
optimizers are `COMPUTATIONAL-EVIDENCE`.  Paper Table-1 values are
`REPORTED`.  A count is called machine-verified only after the independently
known startup count succeeds and the exact coverage construction is checked.

## 6. Pre-registered Rule-7 scope sentence

This campaign sweeps exactly all unit cubics in the orthonormal-probabilist
Hermite span `psi_0,...,psi_3` through their complete even-coefficient disk
`(a0,a2) in [-1,1]^2`, with all compatible odd coefficients covered by the
radius `sqrt(1-dist(0,B)^2)`, over the six full c-bands `[1.0,1.3]`,
`[1.30,1.45]`, `[1.45,1.75]`, `[1.75,3.5]`, `[3.50,4.083]`, and
`[4.083,6.0]`, using ratio-1.02 c-tiles, the stated box floor/depth, panel
ladders through 32768, 256-bit Arb arithmetic, `[0,8]` panels, and the closed
`[8,infinity)` tail majorant; it additionally sweeps only the named
`1e-4`-radius anchor/razor coefficient cells at c-pairs `(1.30,1.45)` and
`(6.0,6.05)`, and it does **not** search c in `[0.993405,1.0)`, the remainder
of `[6.0,12.0]`, c above 12, non-unit cubics, degrees above three, the paper's
C1 or C3 certificates, the splice certificate, Gate-C septic/Krivine scheme
families, or any construction outside the paper's unit-cubic threshold
family.
