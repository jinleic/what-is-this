# REFINEMENT LOG — pre-compute corrections to the frozen pre-statement
# (rule 5/16: corrections recorded beside, original byte-frozen; committed
# BEFORE the first certification compute; the pre-statement commit is 739e63f)

All three points below were found while building the instrument code, before
any certification run. No ladder rung, box, seed, slice rule, cap, or
selection rule is changed. The margin ALGEBRA is completed (point 1); the
P-pos construction is disambiguated (point 2); the quarantine criterion is
restated on sound grounds (point 3).

## R1. The exact remainder bound W_e (completes pre-statement sec 2.2)

The pre-statement's per-equation remainder formula Q_e = sum_r |a0[p,r]|
rho^2 is INCOMPLETE, and its justification sentence "the CP Jacobian is
constant" is WRONG for the unsymmetrized CP form: J entries are the degree-2
products b0[b,r] c0[c,r] etc., which depend on x0. What IS true: the Taylor
expansion of the trilinear monoid g[p,b,c] = sum_r a[p,r] b[b,r] c[c,r]
about x0 TERMINATES — with displacement (da, db, dc):

  g(x0+d) - g(x0) = sum_r [ da bc + a db c + ab dc            (linear)
                          + da db c + da b dc + a db dc       (quadratic)
                          + da db dc ]                        (cubic)

The linear pieces are EXACTLY J(x0) d (the Jacobian evaluated at x0), so the
gate-A bounce c0 = Y g0 with Y = J_S(x0)^{-1} cancels them exactly — the
cancellation argument in the pre-statement stands. The remainder bound per
equation e = (p,b,c) is the EXACT absolute-value sum over the quadratic and
cubic pieces, with per-coordinate displacement radii
  Pa[p,r] = rho if a-coord (p,r) in S else 0, and likewise Pb, Pc:

  W_e(rho) = sum_r [ Pa*(Pb*|c0[c,r]| + Pc*|b0[b,r]|)      (quadratic)
                    + Pb*Pc*|a0[p,r]| + Pa*Pb*Pc ]          (cubic)

(all factors at their exact seed values; a0/b0/c0 are the dyadic seed
entries). R2_c(rho) = sum_e |Y[c,e]| * W_e(rho) exactly as pre-registered.
The containment tests are UNCHANGED in form:
  containment: |c0_c| + R2_c(rho) < rho  for all c (strict; fmpq exact;
               arb-outward re-verified);
  exclusion:   exists c: |c0_c| > rho + R2_c(rho)  (exact fmpq; certifies NO
               root of the square slice system in the named box).
Rationale of the exclusion test: any root x* in the box satisfies
0 = g0 + J_S d + w(d), d = x*_S - x0, |w_e(d)| <= W_e(rho); multiplying by
Y gives |c0_c| <= |d_c| + |(Y w)_c| <= rho + R2_c(rho); the negation
certifies emptiness (this is the standard Krawczyk existence/exclusion pair).
The cubic term is dominated by the quadratic terms at every rung (Pc <= rho
<= 1e-2 << |b0|, |c0| ~ O(1)) but is carried exactly anyway.

## R2. P-pos construction disambiguated

The plant root must lie INSIDE the certified box: the pre-statement's
"+1e-6 * default_rng(seed=20260902) normal" perturbation applied to ALL 247
coordinates would displace F-coordinates too, and F is frozen at seed values
— the exact root would then sit OFF the slice and containment could never
fire (an invalid control, not an instrument failure). Disambiguation: the
plant seed keeps F-coordinates EXACTLY at the drawn factor values (frozen =
root values, zero displacement) and perturbs ONLY the 192 S-coordinates by
1e-6 * default_rng(seed=20260902).normal (S determined by the same QRCP rule
at the plant seed). The box then provably contains the root iff the S-radius
covers max_i |dS_i| <= rho, i.e. rho = 1e-2 on the ladder top covers it with
~3 orders of margin.

## R3. P-pos perturbation radius 1e-6 -> 1e-8 (same seeds, same distribution)

Containment needs rho > |c0|_max ~ ||Y|| * ||g0|| while the quadratic margin
needs R2_c(rho) < rho - |c0_c| with R2 ~ ||Y|| * W(rho) ~ ||Y|| * O(rho^2).
The feasible rung window is roughly
  ||Y||*||g0|| << rho << sigma_min(J_S) / (sum of |factor| scales),
a NONEMPTY but NARROW interval. With a 1e-6 plant displacement, ||g0|| ~
4e-5 and the window closes at the wrong end for the top rungs (|c0| ~ 8e-4
exceeds rho=1e-4, while R2(1e-3) ~ 8e-4 pushes rho=1e-3 out). With 1e-8
(altogether still 8 orders above the float64 arithmetic floor and 6 orders
above the certificate-relevant entry scale), ||g0|| ~ 4e-7, |c0| ~ 8e-6, and
the window rho in [1e-5, 1e-3] is comfortably open on the fixed ladder. The
seed RNGs and distribution are unchanged; only the scale index changes. This
is a control-construction fix, registered before any run; it does not touch
the main system, its seeds, or its ladder.

## R4. Sound quarantine criterion (replaces the pre-statement's P-wrong
## quarantine clause; C-tau-6 becomes the load-bearing counterfactual)

The pre-statement said containment firing on P-wrong (corrupted target)
quarantines the instrument. That is UNSOUND as stated: a containment
certificate on the corrupted system T_cor would be a TRUE statement (an
exact rank-13 CP decomposition of T_cor exists in the box); a true
certificate is never an instrument failure, and T_cor (differing from T_F in
one entry among 192, with 55 gauge-dimension directions frozen) may genuinely
possess nearby roots. The SOUND falsification criteria, registered now:

  Q1 (false-certificate test, decisive): containment firing at ANY rung on
  the C-tau-6 system (tau at r=6) => INVALID INSTRUMENT. Grounds: tau has
  EXACT rank 7 (machine-verified both directions in-repo), so the r=6 CP
  system has NO root anywhere; a containment certificate would be a false
  existence proof. Expected: no containment at any rung AND certified local
  exclusion (disjointness) at the small rungs: at rho -> 0, R2 -> 0 like
  rho^2 while |c0| stays bounded away from 0 (the polished r=6 residual is
  bounded below by the exactly-certified rank gap), so |c0_c| > rho + R2_c
  provably fires for some c at sufficiently small rho on the fixed ladder.
  C-tau-6 is therefore the assignment's "deliberately wrong candidate the
  Krawczyk test MUST reject" — with the rejection being a MACHINE-VERIFIED
  no-root-in-named-box certificate, not merely an absence.

  Q2 (determinism): any (system, box) evaluated twice yields identical exact
  outputs (trivially expected in exact arithmetic; checked once).

  Q3 (arb/fmpq consistency): every fmpq-strict margin comparison must agree
  in direction with its outward-arb re-evaluation; a disagreement means an
  implementation bug => quarantine.

P-wrong runs as registered and its outcome is REPORTED with its own label;
if containment fires on T_cor, the certificate is labeled as being about
T_cor ONLY (it would be an exact fact about a corrupted table, not about
T_F, and it triggers no quarantine and no claim about Route F).

## R5. Main-system pass-window honesty note (no rule change)

With the completed W_e, the main-system PASS window at polish residual rel is
approximately rel < sigma_min(J_S)^2 / (127 * scale^2) (bounce requirement
rho > ||Y||*||g0|| versus margin requirement rho < sigma_min/(sum |factors|)).
At the frozen CP seed the lexicographic sigma_min was 6.1e-5 (float recon);
QR-pivoted selection should be no worse. If the polish stalls at rel ~
1e-3 (as the pre-registration recon found from the raw S_CP seed), NO rung
can pass and the outcome is FAILURE TO CERTIFY — reported with the full
bounce/margin/exclusion data per rung; the adjudication vocabulary is
unchanged.

## R6. Correction entered at first dry-run (BEFORE any certification
## compute; original pre-statement text stays byte-frozen)

The pre-statement's C-tau-7 count "64 equations (4*4*3)" is WRONG: tau is a
3 x 4 x 4 tensor, so the equation count is 3*4*4 = 48 (the "64" was a
4*4*4-style miscount). The correct counts, now verified in-code:
  C-tau-7: 48 equations, 77 unknowns, frozen 29 (48 = 77 - 29; balanced).
  C-tau-6: 48 equations, 66 unknowns, frozen 18 (48 = 66 - 18; balanced).
The MAIN (T_F, r=13) counts 192 = 247 - 55 are CORRECT as registered
(3*8*8 = 192; (3+8+8)*13 = 247).
Also per R3 the P-pos plant perturbation scale is 1e-8 (not the
pre-statement body's 1e-6); R3 was registered before the dry-run, resides in
this log, and the instrument implements 1e-8. Where body and R3/R6 disagree,
R3/R6 governs.
Dry-run facts recorded (all before ladder compute): 3-way table anchor PASSES;
exact residual at the frozen tau_r7 certificate factors = 1.8789e-15 max
|g0|; exact rank of the tau-7 Jacobian at that seed = 48/48 (full row rank),
float sv in [0.29, 8.0e0]-ish; the Krawczyk slice system is well-posed.

## R7. C-tau-6 polish budget note (no rule change)

The 64x66 polish numbers in R6's corrected counts mean the C-tau-6 pipeline
is a 48x66 system; max_nfev/round and the 5-round cap are unchanged.
bounce/margin/exclusion data per rung. The campaign's adjudication
vocabulary is unchanged.

-- OctRankRouteF, 2026-09-01, before first certification compute.
