# Progress ledger

Newest first. One entry per session. Every claim links to its verification.

---

### Session 2026-09-01/02 — 36 lifecycle campaigns, 17 certified, every artifact owner-verified, 2026-09-02

One entry for the whole autonomous session, newest first within the results
index; each item links to the campaign directory and its checksum ledger via
[`RESULTS.md`](RESULTS.md). Every subagent artifact below was verified by the
owner independently: checksum ledgers replayed, prereg binding checked against
the manifest (and both hashes where a prereg was amended), and the headline
mathematics cross-checked in the owner's own kernel with a from-scratch
instrument. Session-wide integrity pass at close: **36 campaigns, 0 checksum
mismatches**, verdicts 17 FROZEN-CERTIFIED / 7 FROZEN-NEGATIVE /
4 FROZEN-INCONCLUSIVE / 1 CRASHED / 1 REJECTED / 1 REHEARSAL, plus the live
kg band, mceliece `t=48` and delcap `n=13` runs.

**Theorems (sharp hypotheses, closed-form counts, exhaustive verification).**
`rs-pe3d/`: H-MINLINE **falsified** at $q=13,s=(2,2,4)$ (24 = 16 lines + 8
diagonals, explicit witness); **T-DGE** certified (spark of a Kronecker product
is the minimum spark; size-$d$ dependent sets are axis fibers iff at most one
factor has spark 2 — sharper than "$d\ge3$"); at $|S|=d+1$ the fiber picture
**collapses** (156 size-4 circuits, none fibers), **Thm-CRIT** $d_A+d_B\le d+3$,
**Thm-FIB**, the **crossing** theorem with count $d_Ad_BC_AC_B$ and its
converse completed by **H2**, the **cross-ratio** law for the all-distinct
channel and its **PGL(2,p)** closed forms, the **$n$-factor** theorem (at most
two varying coordinates at size $d+1$; all fibers for $d\ge4$), the complete
**size-5** classification (degree-$\le2$ bipartite graphs with minor
conditions), the **two-row completeness** capstone (circuit size
$\le r_Ar_B+1$), the mixed and tied three-row layers, and the size-8 residual
explained as **complete intersections of two $(2,2)$-curves**. Twelve campaigns;
the target rests.

**Frontier moves with replayable certificates.** `mm3/`: the fixed-orientation
ladder is now exact at every rung — `paper55` $\sigma^{0/1/2}=55$,
`sun56`$=56$, `mws59`$=\mathbf{58}$ (an explicit 58-addition circuit, one
below the published 59, re-expanded by the owner over all 81 monomials),
`stapleton60`$=60$, `laderman23`$=62$ from the source-locked 1976 primary
(printed basic form 98). `kg/`: tile 3 of $[4.083,6.0]$ certifies at cap
$2^{19}$ (owner recompute agrees to $10^{-35}$); band $[1.75,3.5]$ measured to
close at $2^{17}$ and its campaign is running. `delcap/`: the $q=3,n=12$ box
complete, 4/4 rows beating both published endpoints. `mceliece/`: the
$t=96$ row at $m=12$ with $\alpha$ measured on all 3488 points; `t=48` running.

**First-class negatives.** `omega/`: $\varepsilon$-freezing and kernel-ray
restriction both falsified with measured factors (355.9× and 223× signal; the
straddle is not a box artifact), one Tier-2 exclusion certified, $\varepsilon_0$
proved bit-constant; the 21-dimensional question stays open and no record is
implied. `oct-rank/`: N1 (0/40 seeds at the exact gate) and N2 (F4 blows a
30-minute cap on a system four times smaller than the target, $8\text{M}\times74\text{M}$
matrix, 22.65 GB) — both routes to $\{13,14\}$ closed with costs; the frontier
is unchanged and the absence sentences stand verbatim.

**Corrections to frozen evidence, owner-adjudicated.** The session-3
`mm3/` `transpose_check.py` was found self-circular (its `chain` dict never
populated; its own stdout contradicts its JSON) — the number it supported is
correct and now independently established, the artifact is void. Two `omega/`
provenance records read misleadingly after prereg amendments; corrected at
target level with the rule that amended preregs carry both hashes, which every
later campaign then honoured.

**Process facts worth keeping.** Pre-commit compute was caught and forced into
REHEARSAL/REJECTED restarts three times (mceliece, kg twice); every
subagent-reported "throttling" this session was either a 600 s / 7200 s
`RLIMIT_CPU` in the eval kernel (exit 152 = SIGXCPU) or nice-15 priority
starvation against a nice-10 C++ enumeration at load ~50 — priority is not a
scientific knob and was disclosed wherever changed. Bytecode was never written
into a frozen directory by the owner; one `__pycache__` created by a sibling at
05:28Z inside campB's code dir is recorded as a hygiene defect, not deleted.

---

### `mceliece/` m=12 REACHED at a real NIST cell; alpha measured, not derived; the section-3.6 hole stays BOUNDED, 2026-09-01

Agent `MceliecelM12`, campaign
`mceliece/campaigns/2026-09-01T03-28-10Z_M12ALPHA_DIVONLY/` (commits `7fb44ae`
pre-statement before compute, `30b6c8f` AMENDMENT_1+2 with calibration, `6ac1e44`
anchors + AMENDMENT_3, `fc88126` plants, `240820d` freeze, `6027427` checksum
ledger 17/17, `bd74601` README append). Owner-verified from `verdict.json`.

**`m=12` was UNREACHED and frozen as such; it is now reached, and at a real
Classic McEliece parameter set rather than a toy.** Instance
`(m, n, t, seed) = (12, 3488, 64, 16384)`, i.e. the `mceliece348864` cell.
Owner cross-checks of the frame: `k = 3488 - 12*64 = 2720` matches the artifact,
and `D = n - 2t - 1 = 3488 - 128 - 1 = 3359` matches both the artifact's `D` and
its `max_deg_f`, so the degree gate closes exactly rather than with slack.

**alpha is MEASURED, not derived — the load-bearing distinction.** ADDENDUM 2
pre-registered two paths: measure alpha, or fall back to a derived
**CITED-DEPENDENCY** value if the measurement exceeded budget. The derive fork was
never reached. `alpha_check_kind` reads *"MEASURED (all n=3488 support points,
unmodified m<=11 instrument ops, exact table arithmetic)"*, with
`alpha_grid_points_checked = 3488`, `alpha_failures = []` and
`alpha_identity_measured_full_grid = true`. What licenses "identical instrument"
is the anchor set: 3/3 byte-equal against the frozen gate-A records at
`(11,2048,48,6211)`, `(10,1024,40,5113)` and `(6,64,3,1387)`, with only
`elapsed_s` excluded (declared) and an additive `guards_build` field disclosed. A
derived alpha presented as measured is precisely the defect this target recorded
once before; it did not recur.

**delta is exact, and the one component I suspected of being sampled is not.**
The chain: exact division `3488/3488`, Lagrange unit checks `3488/3488` with no
failures, `offdiag_probes = 200` with `0` violations, assembly `3/3`, degree gate
`3359 = 3359`. The `200` probes initially read like a sample of roughly `6.08e6`
off-diagonal pairs, which would have undercut the exactness claim — but
`delta_scope` states the off-diagonal vanishing is STRUCTURAL, and that is
correct: for a Lagrange basis, `L_i(a_l) = 0` at `l != i` is definitional, since
the numerator carries the factor `(a_l - a_l)`. So the exhaustive content is
`L_i(a_i) = 1` at all `3488` points plus the assembly and degree gate, and the
probes are a redundant control. The **MACHINE-VERIFIED** label holds.

**Counterfactual plants 3/3 REJECTED**, each by a different mechanism: CF-1
duplicate support (`Pi' = 0` at `[10,11]`); CF-2 `f0 + Pi` triggering a
degree-gate ABORT at `2048 > 1951`, with value-invisibility separately CONFIRMED
`0/22`; CF-3 `row0 -> 3*f0` failing the alpha value check (`c^2 = 5 != 3`).

**Scope, stated exactly and not inflated.** This is a BOUNDED FINITE statement
about ONE instance. `t = 48` and `t = 96` at `m = 12` remain **UNREACHED**; no
`n != 3488`, no other orderings, no `m > 12`, no `m = 12` census rows, and no `G`
population exhausted. It is not an `m=12`-complete result and not a NIST
attack-cost verification, and **Apon's section-3.6 hole remains BOUNDED, not
closed**. Cost measured in-process only per rule 17e: per-point alpha `0.1504`
CPU-s, alpha grid `700.53` CPU-s, build `2175.43`, total `2893.37` of the
registered `10800` budget, wall about 51 minutes.

**Six defects disclosed**, including three amendments made and committed BEFORE
the verdict compute (AMENDMENT_1 recording that the measured projection `524.6`
CPU-s beat the `10800` budget so path-1 applied; AMENDMENT_2 correcting the
agent's own CF-2 expectation; AMENDMENT_3 replacing an inert CF-3 mechanism), a
disclosed `1.33x` under-projection of the alpha grid absorbed by headroom, a
`3600`-s tool cap hit that forced split stages with no partial verdict used, and a
briefly overwritten checksum ledger corrected in-run to a complete 17/17.

Named next action: a second `m=12` row at `t = 96` under a fresh pre-statement,
budgeted from the measured rates (build about `2175` s, alpha about `700` s, delta
about `15` s per row).

### `kg/` the pre-registered subdivision is PROVABLY DEGENERATE; the obstruction identifies the real binding constraint, 2026-09-01

Agent `KgCampB`, campaign `kg/campaigns/20260901T051620Z_KgCampB_subdiv_env/`
(pre-statement `d557e70` before compute, skeleton `01bd720`, freeze `f744dd2`,
checksum replay 11/11 OK). Phase-1 FAIL by the pre-registered threshold; Phase 2
correctly NOT run.

**The refinement kg had carried as its named next job was never going to work,
and now that is a certified fact.** The `mq = t^2` panel subdivision recovers
**exactly zero** on Campaign A's failing tile-3 cell: `V - U = 0 +/- 2.36e-73`.
Mechanism, owner-checked by hand: the crossing set
`{s > 0 : odd2 * (5/2 - s^2 + s^4/6) = c_L^2}` has its unique root above
`sqrt(3)` at `s* = 45.703407790989350730` (certified bracket) — far outside the
`[0,8]` range — and no in-range root below `sqrt(3)` exists because
`(5/2) * odd2 = 5.93e-5` is dwarfed by `c_L^2 = 17.21`. Reconstructing
independently: `odd2 ~ 2.372e-5` forces `(5/2 - s^2 + s^4/6) ~ 7.26e5`, hence
`s ~ 45.7`, matching the certified root. The certified second-hinge kill fires
only on panels where A's bound is already zero. So the mechanism targets
structure that is absent from the domain — a **quantified degeneracy**, not a
tuning shortfall.

**Bit-exact cross-campaign reproduction.** B re-derived A's tile-3 margin as
`-1.76291787028287693244667476546e-5 +/- 4.20e-35`, identical to A's frozen
value, with `recovery_vs_A = -4.19e-35` (i.e. zero) and the reconciliation gate
reporting `reproduces_A_frozen = true`. Two campaigns, two instruments, the same
number to the last digit.

**The obstruction is constructive, and it reframes A's failure.** The frozen
decomposition, with the owner re-checking both additions:
`D1 = -1.76291787028287693e-5` at cap `2^17`; `D2` fold recovery `= 0`;
`D3` panel-sup gain `2^17 -> 2^19 = +2.20922696204705e-5`, and
`-1.76291787028287693e-5 + 2.20922696204705e-5 = +4.4630909176417e-6`, matching
the frozen `D3_margin_cap_2p19 = +4.46309091764170239642742026145e-6` exactly;
a further `2^19 -> 2^21` gain of `+5.52295164708941e-6` gives
`+9.98604256473110786e-6`, again matching. `D4` shows a coefficient-space
disk-feasibility diagnostic with leverage `4.94055685681106e-4` = **28x the
deficit**, correctly **NOT adopted** because it needs a sound joint
`(even, odd)` re-derivation first — the firewall held. `D5` tail
`3.6353172540086755822e-13` is identical across all variants. Measured
single-cell cost at cap `2^19`: about `68` s CPU.
**Consequence: the cell is not intrinsically hard.** Panel-sup gains alone
(`2.76e-5` across `2^17 -> 2^21`) exceed the deficit, the frontier leaf is depth
`14/40` and not floor-bound, and A's stop reason was the panel cap. A's own
disclosed addendum had weakened the cap from its originally pre-registered
`2^19 = 524288` down to `2^17 = 131072` — a weakening in the SAFE direction,
since a weaker instrument can only fail to certify and never manufacture a PASS,
but it is what cost the band.

Standing status unchanged by this campaign: band `[4.083, 6.0]` remains **FAILURE
TO CERTIFY** at tile 3 of 49, A's tiles 1-2 PASSes stand, the `kgcloseB3`
retraction **STANDS**, and nothing here bears on the paper. One mid-campaign
defect disclosed (a first-hinge `W`-drop) was caught by the agent's own
reconciliation gate and fixed before any verdict existed.

Next: Campaign C dispatched with the MINIMAL attribution-clean change — tile 3
re-run entirely at cap `2^19` with every other setting byte-identical and one
knob moved, gated on a bit-exact reproduction of A's cap-`2^17` margin first,
then tiles 4-49 with a pre-registered global budget, incremental per-tile freezes
and tile-granular resumability. B's `D4` diagnostic may be adopted only after a
sound joint `(even, odd)` derivation that admits the paper's kernel point and
passes the same counterfactual battery; a 28x-leverage diagnostic is not a bound
until it survives that.

### `kg/` band [4.083, 6.0] under the corrected envelope: FAILURE TO CERTIFY at tile 3 of 49; the retraction stands, 2026-09-01

Agent `KgEvenBoxRerun`, campaign
`kg/campaigns/20260831T231500Z_KgEvenBoxRerun_campA_corr_env/` (commit chain
`e4e8416` -> `4f4594d` -> `c72be46` -> `f6c8c36` -> `d3a1a44` -> `0f494f3` ->
`8d2cb26`). Every number below was owner-re-derived from
`logs/startup_controls.json` and the per-tile artifacts.

**This is the outcome a SOUND envelope should produce where an unsound one had
manufactured a PASS.** The predecessor `kgcloseB3` band PASSes were retracted
because its `j_integrand_bound` capped `|e| <= 1` while the paper's own D.4
reproducing kernel reaches `|e(0)| = sqrt(3/2)`, and because it used 3-D
circumradius inflation instead of the paper's 2-D even-box form. Both defects are
corrected here, and the corrected instrument does NOT reproduce the retracted
PASS. The retraction therefore **STANDS** (rule 5); no new PASS has been earned.

**Verdict: tiles 1-2 PASS, tile 3 OPEN, band FAILURE TO CERTIFY.** Tile 1
`[4.0830000000, 4.1156640000]` margin `+1.04499016763434518123845756638e-6 +/-
4.27e-36`; tile 2 `[4.1156640000, 4.1485893120]` margin
`+4.71904853878616890385690164962e-6 +/- 3.07e-36` — both strictly positive with
radii near `1e-36` (**MACHINE-VERIFIED**, Arb 256). Tile 3
`[4.1485893120, 4.1817780265]` has an open leaf at margin
`-1.76291787028287693244667476546e-5 +/- 4.20e-35` with `95/251` leaves certified
at depth 14 and the panel cap `131072` reached; frontier
`a0 in [-0.830729166666666666666667, -0.829427083333333333333333]`,
`a2 in [+0.558593750000000000000000, +0.559895833333333333333333]`. Per the
pre-registered stop-at-first-failing-tile rule, tiles 4-49 were NOT run. The cell
is a disk-ring where `r_o -> 0` and the binding cost is panel-sup `|e|`-hinge
slack. Correctly scoped by the agent and confirmed by the owner: the negative
margin is a failure of the panel-sup envelope at the caps — it is **not** a lower
bound on `J - d` and **not** a paper refutation.

**The startup battery is the strongest instrument validation in this target, and
it validates the verdict predicate itself, not just the envelope.** 16/16 PASS:
`e(0) = 1.2247448713915890490986420373529457 +/- 4.02e-36`, equal to `sqrt(3/2)`
to every shown digit; the envelope ADMITS that kernel point and is TIGHT there
(slack `1.87e-96 < 2^-260`); the cap-at-1 counterfactual is rejected with the
clipping quantified as `sqrt(3/2) - 1 = 0.224744871391589049098642037353`; the
3-D-inflation counterfactual confirms the OLD form was looser
(`1.65423992279850488104968` versus corner-exact `1.32426406871192851464051`,
slack `0.32997585`); the 2-D envelope is corner-exact at `s = 0`; the authors'
`d3h` certificate is byte-verified (`sha256 80e945589b...`) and its `B3` is
imported UPWARD (`14.4424366466397645700000000000` above the certificate's lower
end); and the margin predicate is validated three ways — a straddling ball is
rejected, a separated ball accepted, an upper-above-`d` ball rejected. The cited
tail `4.5756855097934878218797010998724522e-6 +/- 4.62e-41` is consumed as a
cited constant, independent of the new envelope.

**One precision correction to the agent's own harness note, recorded rather than
adopted.** The report attributes six earlier detached launches producing zero
verdicts to "an env artifact ~40-50x". That cause was never established: per rule
17e the `ps`-derived rate readings on this box are unreliable, no `/usr/bin/sample`
stack was retained before those runs were abandoned, and the in-session tile walls
(`2665`/`3455`/`863` s) mean a detached launch would need over 44 minutes to
produce even tile 1, so "zero verdicts" is consistent with several causes. The
honest record is: six detached launches yielded nothing, the cause is
**UNDIAGNOSED**, and the verdict run completed in-session. Campaign B's dispatch
carries an explicit instruction to retain a stack before killing anything.

Also disclosed and accepted: pre-statement addenda superseding ratio `1.02` with
`1.008` (the former is arithmetically impossible — certified `-1.06e-3` at the
kernel) and panel cap `524288` with `131072`, both before the verdict run.

Next: Campaign B is dispatched — the pre-registered `mq = t^2` panel subdivision,
validated FIRST on exactly this tile-3 frontier (mechanism validation on the
known-failing cell, with the deficit accounting reported whatever happens), and
only on a strictly positive margin does it continue tiles 3-49 of the band. The
three never-run bands `[1.0, 1.3]`, `[1.45, 1.75]`, `[1.75, 3.5]` are Campaign C,
deliberately withheld: deploying an unvalidated mechanism on fresh bands before
testing it where the current instrument demonstrably fails is backwards.

### `oct-rank/` Route F: FAILURE TO CERTIFY with a certified local negative; the instrument proved it can discriminate, 2026-09-01

Agent `OctRankRouteF`, campaign
`oct-rank/campaigns/2026-09-01T04:30:00Z_routeF_kraw/` (git `739e63f`
pre-registration before compute, `15e1c95`, `77fbe24`, `3074618` freeze,
`1245009` hygiene, `d505584` README append).

**The agent caught a rule-14 defect in the plan the owner relayed, before
spending compute.** The inherited `af-triple13-existence` sketch specified a
"20-variable symmetric-real / symmetrized-Jacobian" Krawczyk. That instrument
cannot establish the claim: a gradient-system Krawczyk certifies CRITICAL POINTS
of the residual, a positive-residual local minimum is not a decomposition, and no
interval method can prove a residual is exactly zero. The agent documented the
diagnosis in the pre-registration and substituted an existence-carrying
instrument instead of proceeding. Owner authorized the swap with four conditions:
state the square-system counts, plant a genuine NEAR-MISS rather than an
obviously-bad candidate, run `tau` in BOTH directions including required failure
on the known-impossible case, and fix the outcome vocabulary in advance. All four
were met.

**Counts balance, owner-checked:** `247` parameters (`13` rank-one terms
`x 19`) minus `55` frozen gauge coordinates `= 192` equations `= 8 x 8 x 3`
tensor entries. A square system with an isolated root is what existence requires.

**Result: FAILURE TO CERTIFY, with a certified local negative attached.** The
pre-registered `5 x 5000`-nfev polish bottoms at relative `8.816619e-06` —
roughly four orders above the certified window, which is empty at that residual
because containment needs `||Y g0|| < rho` while the quadratic margin needs
`rho <~ 1.9e-5`-ish. Containment is `0/192` at all 15 fixed rungs
(`1e-2 .. 1e-12`), and certified exclusion grows `44/192` at `rho=1e-4`,
`190/192` at `1e-5`, and COMPLETE `192/192` for every `rho <= 1e-6`. So the
campaign certifies that no exact on-slice rank-13 decomposition lies within
`1e-6` of either polished candidate (**MACHINE-VERIFIED**, exact `fmpq` with
outward-`arb` re-verification of every passing comparison). Correctly scoped: this
is **not** a rank-13 impossibility — off-slice roots are not excluded, exact
transversality is `rank(J_F) = 55/55`, and nothing global is claimed. At the
partial `1e-5` rung the two non-excluded coordinates are `c`-block `(7,6)` and
`b`-block `(4,8)`.

**The instrument was proved able to discriminate, which is the whole point of
rule 14.** Owner replayed `rf_replay.py` from the frozen bytes: checksums OK,
15/15 PASS, 42.6 s CPU. `tau`-7 containment climbs `1/48` at `rho=1e-1` to
`45/48` at `1e-2` to `48/48` at `1e-3` — it ACCEPTS a known-true case and thereby
independently re-certifies the frozen `tau_r7` certificate as an existence
statement. The globally infeasible `tau`-6 plant is REJECTED with complete
certified exclusion at `rho <= 1e-3`; the synthetic rank-13 positive control is
ACCEPTED at `1e-4`; and the near-miss corrupted-target candidate (entry
`[0,0,0]=2`) is REJECTED everywhere. Five cases, five correct directions:
the instrument does not certify proximity into existence.

**Verdict and absences, stated with their exact meaning.**
`rank((L_1, L_i, L_j))` stays **OPEN in {13,14}**, the same certified interval as
before the campaign. Absence of a 13-witness is NOT evidence for 14; absence of an
impossibility argument is NOT evidence for 13; and the certified local exclusions
are not evidence for 14 either. The published window
`18 <= R_R(T_O) <= 25` is untouched and no published work is refuted.

**Four defects disclosed, every one caught by a pre-registered control before
certification data existed:** a run-1 target-convention error (raw slices versus
blockdiag-conjugated target) caught by the seed-reproduction anchors, with seeds
then reproducing bit-exactly on the corrected target; an entrywise-false first
anchor-block formulation caught by the instrument's own startup halts; a singular
`tau`-7 natural slice re-run under the registered `qrcp` rule; and two count/scale
slips (`64` vs `48` equations, P-pos `1e-6` vs `1e-8`) corrected pre-compute in
`REFINEMENT_LOG` R3/R6/R7. Cost reported via `time.process_time()` only, per rule
17e — `415.5` s instrument, `42.8` s replay, no `ps`-derived numbers.

Named next actions, both requiring owner re-scoping per rule 16: N1 exact
CP-completion downward, N2 Groebner infeasibility upward.

### `rs-pe3d/` Theorem T re-derived in corrected form with a commit-anchored pre-statement; owner verification by a second instrument, plus a sharper new finding, 2026-09-01

Agent `RsPe3dMW` (two campaigns, then a clean budget handback). Campaign M:
`campaigns/2026-09-01T02-33-35Z_mech2_mechanismCorrected/` (git `29b0dfe`
pre-statement committed BEFORE compute — this closes the predecessor campaign's
**REPORTED**-grade chronology gap — then `e117ae2`, `d45174a`). Campaign W:
`campaigns/2026-09-01T20-30-00Z_win4_pn4/`, pre-statement `8114a8f`, frozen
mid-flight at budget with validation and resume complete.

**Theorem T, corrected form (HUMAN-AUDITED derivation + MACHINE-VERIFIED
controls, 16/16 exact `F_q`).** (1) For EVERY support with `|S| < d`,
`V cap F^S = {0}`. (2) There EXISTS a support of size `d` carrying a codeword.
(3) NO per-support claim at `|S| >= d` — clause C-1, honored explicitly. Census
semantics: the union window is empty iff `w < d`, first nonempty exactly at `d`.
At `t=(1,1,1)`, `d = min_i s_i`, so H-GATE's min-form is a **theorem** on that
slice, and the `(4,5,7)` class the P2 battery could not test is inside its grip.

**Owner audit of the mathematics, by hand.** L5 is the load-bearing step and it
is valid: build `y = phi_0 (x) phi_1 (x) phi_2` with each `phi_i` killing `C_i`,
so `y` annihilates every generator of `V` (each has one dead factor); the
existence of each `phi_i` needs exactly
`pi_i(e_{x_i}) not in span{pi_i(e_a) : a in T_i \ {x_i}}`, which L3 supplies
because `|T_i| <= 1 + (w-1) = w < d <= d_i`; then
`0 = <y,v> = v_x != 0`. Induction-free and correct. Restated, Theorem T says
`d(V) = min_i d(C_i)`: the `<=` direction is trivial from each slice, and L5 is
the real content — no cancellation across the three summands lightens a word.

**Owner verification by a DIFFERENT instrument.** The agent builds `V` as an rref
of stacked line classes; the owner built it as
`V = ker(H_0 (x) H_1 (x) H_2)` from GRS parity checks and tested emptiness via
`dim(V cap F^S) = |S| - rank(pi[:,S])`. Reproduced exactly: PN4 `dim V = 68`; the
`w <= 3` census `457,450` supports with `0` nonzero (agent: 862 s, owner: 9.4 s);
`35/35` direction-0 minimum lines carrying codewords. Counting cross-checks all
consistent by inspection: `457,450 = 140 + 9,730 + 447,580`, `43,744 = 64 + 2,016
+ 41,664`, `85,400 = 80 + 3,160 + 82,160`, `C(61,2) = 1,830`, `C(63,2) = 1,953`,
corner censuses `14` and `41`. The owner also EXTENDED the anchoring: every
emptiness control in the campaign was `t=(1,1,1)`, so parts 1 and 2 were tested
at `(13,(5,4,4),(2,2,2))`, `(13,(4,4,4),(2,1,3))` and `(31,(5,5,6),(3,2,2))` with
random `Lambda` — all PASS, so the theorem's machine anchoring is no longer
`t=1`-only. Two false FAILs appeared first and were the OWNER's bug (a leaked
loop variable made `d_i` wrong whenever the `t_i` differ), the same class as three
of the agent's five disclosed watch-points; the dependent triples my buggy run
flagged were exactly part 2's predicted witnesses.

**New owner finding, sharper than clause (3): "H-MINLINE".** At weight exactly
`d`, exhaustive exact censuses return candidate sets SET-EQUAL to the minimum
lines in argmin-`d` directions, with ZERO non-line dependencies — PN4 `35 = 35`
over all `C(140,4) = 15,329,615` supports; `(17,(4,4,4))` `48 = 48` with all three
directions minimum; `(41,(4,4,5))` `40 = 40` with directions 0 and 1 only and
direction 2 correctly absent because `s_2 = 5 > d`. Status stated exactly:
**MACHINE-VERIFIED** as a finite statement about those three instances,
**CONJECTURE** as a general claim. If proved it upgrades clause (2) from
existence to exact identification. The proof attack was handed to the target with
its binding step named — run L5 in reverse, where `|S| = d` forces every valid
axis-assignment to leave some axis with `|T_i| >= d_i`, which pins all of
`S \ {x}` onto one argmin axis with `d` distinct coordinates; whether that forces
a full line is where it closes or fails. It is pre-registered as its OWN campaign,
deliberately not folded into W.

**Instrument evidence, and a methodological win.** Two independently built kernels
agree `15/15` on the 13 frozen gate-B candidate tallies plus `PN6 3486` and
`PP2 917`, with B11's `2,304,200`-support census landing in 75 s against the
frozen 3,449 s. The owner ruled two process questions during the run: the census
and the pre-registered ladder are independent computations, so they run
CONCURRENTLY rather than trading off (the agent had framed it as a 40-minute
delay); and the pre-registered V-basis ladder was NOT killed when an interim
result made its remaining rows look redundant, because **pre-registered evidence
must not be dropped on a prediction of its own outcome** — self-sealing reasoning
that costs nothing if right and destroys the only informative evidence if wrong.

Next: a finishing agent holds W's resume protocol (parity primary over the
remaining supports, element-wise candidate-SET comparison against the V-basis
full-space census, exact `rho` closure over the candidates, P and H-GATE verdicts,
ladder rows folded into the freeze). One correction carried into that dispatch:
the first landed line-support hit is an instance of Theorem T part **2**
(existence at weight `d`), not part 1, which is emptiness strictly below `d`.

### `oct-rank/` Route A CLOSED-NEGATIVE with the Strassen 1983 primary read first-hand; the 16-route dies twice, 2026-09-01

Agent `OctRankRouteA`, campaign
`oct-rank/campaigns/2026-09-01T03:21:00Z_routeA_prov/` (commits `dd40d2b`
pre-statement pre-compute, `076a4dd` freeze, `c8defdc` Route F status note).

**The provenance gate that killed the predecessor is now genuinely cleared.**
Route AF had closed Route A using Landsberg's survey and Koiran's restatement
while flagging Strassen primary provenance **UNESTABLISHED** — and the earlier
`16 >= at n=38` claim had failed precisely on that unread source. This campaign
obtained the primary: ScienceDirect 403 (even via a 2024 Wayback capture), EuDML
searched and recorded as a **content-layer** negative (rule 17d applied by the
agent unprompted — the page text was read, not a status code), then
scholar.archive.org -> CORE 82110360 -> a 2019-04-16 Wayback
`application/pdf` capture. Owner verification of identity, from the retained text
layer rather than the report: title page reads *Rank and Optimal Computation of
Generic Tensors / V. Strassen / Institut f[ue]r Angewandte Mathematik /
Universit[ae]t Z[ue]rich*, running head `52/53:645-685 (1983)`, Elsevier, with a
period-correct bibliography (Alder-Strassen 1981, Atkinson-Lloyd 1980, Kruskal
1977); retained publisher metadata matches DOI
`10.1016/0024-3795(83)80041-X`, first page 645. PDF sha256
`f8481709de3a93e5747ce45e5c19aae81225139f39db347ef373a82df5bd0376`
(1,936,012 bytes); text layer sha256 `b138a8b8...`.

**Verdict: the theorem family caps at 12, below our certified floor of 13.**
Theorem 4.1's display carries Strassen's border-rank `R` with underbar while the
proof's `(4.2)` step establishes the plain-rank form, and the form carries a
**one-half** factor: `rank >= n + (1/2) rank(comm)`. Applied to the 3-slice
`(L_u, L_v, L_w)` families this gives `8 + 8/2 = 12`. Owner strengthening,
communicated to the agent: the closure does not need the campaign's `M^2`
identity at all, because `comm` is an `8x8` matrix so `rank(comm) <= 8`
trivially, capping the family at `12 < 13` unconditionally.

**The naked `n + rank(comm) = 16` reading dies twice.** First on provenance: it
does not exist in the primary text — Section 4 contains only Theorem 4.1 and
Corollary 4.2, Blaeser is uncited (the paper is 1983), "commutator" is never used
as a name, and Lickteig appears only as an independent announcement. Second, and
independently of any source, on arithmetic: applied to the `tau` control it would
assert `4 + 4 = 8` while `tau`'s true rank is `7`, and a lower bound cannot exceed
the true rank. The `tau` control was mandated in the dispatch precisely because it
had already refuted an earlier reading, and it did so again. Corollary 4.2 is
inapplicable to the octonions because its associativity clause fails — octonions
are non-associative, a definitional check.

**Controls, all PASS (47/47, exact `fmpq`, exit 0):** `L_conj(u) L_u = N(u) I`
exactly on all eight basis units and a generic triple; generic Gram determinant
`39 > 0`; `tau` bound `4 + 4/2 = 6 <= 7`; `tau-box-s` bound `8 + 8/2 = 12 <= 13`.

**A label nuance the owner flagged — and then RETRACTED, because the owner was
wrong (recorded inline per rule 5).** C3 verifies
`M^2 = -4 detGram/N(u)^3 I` exactly at five basis triples plus one generic
triple, and the owner objected that six exact points cannot establish a
polynomial identity in the 24 entries of `(u,v,w)`, citing the target's earlier
`173/173`-is-not-a-quantifier defect. The agent pushed back with evidence rather
than complying: the universal statement is proved in a DIFFERENT frozen campaign,
`campaigns/2026-08-31T08:02:18Z_routeAF/commutator_derivation.txt`. The owner then
audited that derivation line by line and it is valid:
`L_x^T = L_conj(x)` and `L_conj(x)L_y + L_conj(y)L_x = 2<x,y>I` are **bilinear**,
so basis-coefficient anchors plus bilinearity do give them for all `x,y`; from
`<u,z> = 0` one gets `L_z L_conj(u) = -L_u L_conj(z)`, hence `P^2 = -alpha I`,
`Q^2 = -beta I`, `PQ + QP = -2 gamma I`; therefore `C = 2(PQ + gamma I)` and
`(PQ + gamma I)^2 = -P^2 Q^2 + gamma^2 I = (gamma^2 - alpha beta) I`, giving
`C^2 = -4(alpha beta - gamma^2) I`; the Gram block form supplies
`alpha beta - gamma^2 = detGram/U^3`, so `D^2 = -4 U detGram I`, whose scalar is
strictly negative exactly on linearly independent real triples — so `rank(D) = 8`
**universally**, as a human proof with bilinear machine anchors. The six C3 points
are cross-check anchors, not the basis of the quantifier. The owner's flag is
withdrawn; the agent additionally restricted C3's own wording to the six swept
triples, which is fine but was not required. Independently, the same derivation
reproduces the `1/2` factor through Landsberg's `rank([...]) <= 2(r - b)` with
`b = 8`, giving `r >= 12`, and the `tau` calibration `b = 4`, comm rank `4`,
`r >= 6` against a true rank of `7`.

**Scope preserved.** S3 stays OPEN at 13 with gap exactly 1. The published window
`18 <= R_R(T_O) <= 25` is untouched, no published work is refuted, and nothing
required escalation. Route F was NOT computed — the agent hit its budget cap right
after the freeze and chose to freeze a status note rather than half-launch a
certification it could not finish, which is the correct call. That note is the
handoff and a fresh agent now holds the `af-triple13-existence` plan: an exact
20-variable symmetric-real Krawczyk existence attempt seeded from the two frozen
numerical candidates, with a `tau` control and a counterfactual plant required.

### `omega/` branch-and-bound feasibility gate: SPREAD, frozen as a quantified obstruction; stage 2 never opened, 2026-09-01

Agent `OmegaPenaltyProbe`, campaign
`omega/campaigns/2026-09-01T03:15:00Z_OmegaGateCBnBGate_9d2e4a7c/` (pre-statement
`f8b9ca7` committed BEFORE compute, GATE-Z amendment `532781d` pre-verdict, freeze
`4ba4c38`, manifest `ffae3cf`, README append `90cbe00`). The gate was designed so
the campaign could only spend budget if the width was concentrated; it was not,
so it spent none.

**Why the gate existed.** Closing the root box's `73.85977508485684x`
width/signal by uniform refinement needs about `log2(74) = 6.2` bisections per
binding coordinate, i.e. `~2^(6.2k)` subboxes over a `k`-dimensional binding set
in a 21-dimensional kernel. Concentration of the width was therefore the
load-bearing empirical question, and the threshold — at most 6 named coordinates
carrying at least 90% of the L1 width — was fixed in the pre-statement before any
measurement.

**Measured, and owner-re-derived from `stage1_decomposition.json` rather than the
report.** The per-coordinate widths sum to `113.66562530936517`, matching
`GATE_Z.live_l1` exactly (difference 0). Sorted shares run from `4.08%` down to
`0.58%`: the top six carry only `21.8493%`, `40` of 45 coordinates each carry at
least `1%`, `23` carry at least `2%`, and **`k* = 36` coordinates are needed to
reach 90%** — reproduced independently. Threshold FAILS decisively (36 > 6), so
the pre-registered branch is SPREAD. Implied cost at that binding-set size:
`2^(6.2 x 36) = 2^223.2 ~ 1.548e67` subboxes (`~3.9e33` on the agent's milder
reading). Frozen verdict sentence: branch-and-bound over D at this enclosure
technology is **INFEASIBLE AT THIS COST**. Stage 2 never opened and zero subbox
compute was performed — the threshold was honored rather than renegotiated after
seeing the numbers (rule 16).

**Straddle quantified.** Every R branch's minimizer identity is undecided (three
candidate slots live on 6 of 7 branches); the widest candidate interval span is
`R_glob[0]` at `3.2123059392596964e-05`, which owner arithmetic confirms is
`203.1x` the `1.5816497000997742e-7` signal — consistent with the frozen
CandidateWitness picture. One prose correction: the agent described the widths as
"nearly uniform"; they actually span `6.99x` (`0.6639861348561744` to
`4.640574461822501`). The verdict is unaffected, because what defeats
branch-and-bound is diffuseness, not uniformity.

**Disclosed amendment, accepted.** GATE-Z was pre-registered as bit-equality but
attempt 1 measured `1.2e-14` last-digit replay drift, making bit-equality
unachievable; the tolerance was amended to `1e-12`/`1e-6` in a pre-verdict commit
with the live arrays as single source of truth. The achieved agreement is
`3.552713678800501e-13`, inside the amended tolerance by a factor of ~3, with
gradient endpoint deviations `1.2e-14` and `1.78e-15`. The amendment was declared
before the verdict and changed no decision threshold.

**Scope of the negative, stated exactly.** The obstruction binds *this* enclosure
technology over `D = {delta : A.delta = 0 exactly, |delta|_inf <= 1e-7}` on the
region-0 glob dist block. It does **not** bound a different certified-gradient
technology, nor analytic branch-pinning — those remain the honest openings. Gate
C stays **PARTIAL**; the 21-dimensional kernel question is **UNCHANGED and
OPEN**; the `12.808x` guardrail is carried verbatim, so nothing here bears on the
published `2.37155181` or the record `2.371177`.

### Phantom 40x CPU clamp: three agents and the owner chased a broken measuring instrument, 2026-09-01

No research content changed; the cost was hours of agent time and one owner
diagnostic prescription that was itself wrong, so it is recorded as a
methodological failure rather than dropped.

`KgEvenBoxRerun` and `RsPe3dMW` independently reported their compute pinned at
`~2.4%` CPU, cross-confirmed each other, and began restructuring around it - kg
moving to foreground chunked execution after three zero-verdict launches, rs-pe3d
killing and relaunching a ladder. Owner tests appeared to show no clamp
(foreground `1.07` cores, async `0.999`, async while printing 20,000 progress
lines `0.970`, file-redirect `0.912`, plus delcap's historical 19-minute async
run at `0.9994`), which refuted the two cheap hypotheses - background QoS via
`taskpolicy -b`, and stdout pipe backpressure, the owner's own leading candidate.
The owner then prescribed `delta_cpu/delta_wall` from `ps -o time` as the
discriminator, and both agents dutifully measured `0.044` and `0.022` cores with
it.

`RsPe3dMW` settled it with a stack: `/usr/bin/sample` on pid 85270 showed
`2287/2287 ms` inside numpy `LONG_remainder` under its rref loop - a full core,
burning continuously - while `ps` time deltas read `0->1 s` across repeated 25 s
windows. Throughput agreed with the sample and not with `ps`: 9/13 exact gate-B
rows in ~65 s wall is unreachable at 0.04 cores. **There was never a clamp.** The
owner's prescribed instrument was the defective one; the owner's own A/B tests
were right only because they used in-process `resource.getrusage`, never
touching the broken path. The tell was visible early and missed: four unrelated
processes reporting *exactly* `2.4%` meant the instrument was constant, not the
load.

Two artifacts survive. Rule 17e (`README.md`): a rate is a claim about an
accounting source, so name it - in-process `ru_utime`/`process_time` first, then
`/usr/bin/sample` thread time, then throughput against calibration; never `ps`
deltas for framework-Python children; and never restructure a campaign around an
unexplained performance observation. Second: the causes of both original stalls
are now **unrecoverable**, because each agent killed its process before taking a
stack - the evidence died twice in one session, which is why the rule mandates
retaining `sample` output before any kill.

### `omega/` gate C V1 CLOSED under the corrected v11 interval core; LP stage returns FAILURE TO CERTIFY; my dispatch premise was stale, 2026-09-01

Agent `OmegaPenaltyProbe`, campaign
`omega/campaigns/2026-08-31T02:40:00Z_OmegaGateCFinish_v1_4c8a2f1b/` (commits
`ec16efe` pre-statement before compute, `d72d8d6`, `c807f74`, `23c1451`,
provenance addendum `13c27d8`). Every number below was owner-re-derived before
this entry.

**V1 PASSES.** Slope replay raw `2.3715538358350807` vs frozen raw
`2.3715538358350803`: residual `4.440892098500626e-16`, which `math.ulp`
confirms is **exactly one ULP** at that magnitude, against a `5e-13` gate. All
seven R branches, `R_sum=2.8170035674609757` and `M_low=2.0942543887102634`
are float-equal to the frozen record (**MACHINE-VERIFIED**). Rule 14 checked at
source, not from the report: in `run_v1_v11.py` step 1 `importlib`-loads the
frozen `stage_b_rung2.py` from the 2026-08-30 campaign while step 2 calls
`gate_c_slope_pass.run_slope_pass`, with SHA-256 asserts pinning the slope-pass
and endgame sources — the two paths are genuinely distinct, so their agreement
is informative rather than tautological.

**The penalty diagnostic I ordered is answered, and the answer is (b).** The
Lemma-1-off control shifts the raw by `2.029648699330977e-6`, reproduced here to
the last digit; no `~2.9e5` amplification exists anywhere in the chain, so the
`lam_sum` penalty path never carried a 0.9-scale residual.

**LP/ANY-y stage: FAILURE TO CERTIFY**, correctly labelled. The rigorous Arb
enclosure `[-5.841014555620028e-6, +5.841014555620028e-6]` straddles zero at
`73.85977508485684x` the `1.5816497000997742e-7` signal; the enclosure is
**MACHINE-VERIFIED**, the ratio is **INFERENCE**, and the midpoint LP candidate
`-1.5773329015168567e-7` stays **COMPUTATIONAL-EVIDENCE**. Owner arithmetic
checks: `width = 2 x swing` exactly, and the swing is outward of
`L1 x radius` by `5.5e-20` — the conservative direction. Gate C stays
**PARTIAL**; the 21-dimensional kernel question is **UNCHANGED and OPEN**; six
LP probes cannot certify 21 dimensions. Guardrail re-verified rather than
quoted: the gap to the published `2.37155181` is `2.0258350805768544e-6`, i.e.
`12.808x` the best available signal, and even the full `|LP min|` would leave
`1.868e-6` above published — this route can never produce a record.

**My dispatch premise was stale, and that is my error.** The dispatch ordered a
diagnostic against `Om_raw = 3.272425321778391`, residual
`+0.9008714859433109`. That premise had already been retracted at
`omega/README.md:590-599` BEFORE the dispatch went out: the advertised source
(SHA-256 `4504365a…`) raises `NameError: name 's_const' is not defined` at line
191 before its first aggregate checkpoint, so those values have no reproducible
producer and survive as **REPORTED** history only. I carried a dead number into
an instruction as live. The agent detected it, disclosed it instead of quietly
skipping, and executed the handoff order on the live tree anyway — which closed
the one genuinely open item: no slope-replay V1 had ever run under the corrected
v11 core. This is the fourth owner-premise defect today, after the `rs-pe3d`
`min(s)<=3`/PN4 conflation, the `kg` `paper_full.txt` path, and my own reading
of an ECCC soft-404 as a status code. The pattern is unchanged and now
well-evidenced: the least reliable part of an owner instruction is its stated
reason, and agents are right to treat every premise as a hypothesis.

Bookkeeping defects, both disclosed by the agent and now recorded append-only in
`omega/README.md` (`13c27d8`, prefix byte-identical, frozen dir and its
`checksums.sha256` untouched, directory NOT renamed): the campaign directory
name carries UTC `2026-08-31T02:40:00Z` while the run executed
`2026-09-01T02:47-02:53Z` (~24 h stale, ordering vs the 08-31T07:57 GateCDiag
campaign now spelled out), and aborted attempts' stderr was overwritten by
successive redirects — forward rule is per-attempt stderr under distinct
filenames.

Next named action, unchanged and requiring its own pre-registered campaign:
certified subdivision (branch-and-bound) of the box so each subbox pins its
R-branch and Lemma-1 structure, or a tighter certified gradient enclosure. The
`12.8x` guardrail is carried verbatim into that campaign.

### Monday 2026-08-31 frontier window swept and owner-re-derived; zero collisions; third coverage defect found in our own instrument

Window closed. Agent `SweepAug31` swept it; every load-bearing number below was
re-derived by the owner from the live sources before this entry, not taken from
the report. Raw non-authoritative notes:
[`docs/scan-raw/frontier_2026-08-31_arxiv-cscc.md`](docs/scan-raw/frontier_2026-08-31_arxiv-cscc.md),
[`…_arxiv-csit.md`](docs/scan-raw/frontier_2026-08-31_arxiv-csit.md),
[`…_eccc.md`](docs/scan-raw/frontier_2026-08-31_eccc.md).

**arXiv** (Atom `published`, the v1 submission timestamp, bounds the sweep — not
announcement day, not `updated`): `cs.CC` totalResults 16, max published
`2026-08-31T12:04:08Z` (`2608.30664v1`); `cs.IT` totalResults 37, max published
`2026-08-31T11:03:56Z` (`2608.30586v1`). Owner re-query reproduced both counts
and both maxima exactly. **ECCC** (the PUBLICATION date bounds it, cutoff
`2026-08-30T15:05:01Z` = TR26-162's `30 Aug 2026 18:05:01 +0300`; the paper date
differs by a day and report numbers are not chronological): **zero new reports**,
feed `lastBuildDate` = max `pubDate` = TR26-162, no `TR26-163+` present.
**Ceiling after sweep: TR26-162, unchanged.** Triage of all new items against all
seven targets: **no collision, no gate affected, no target opened.** Two genre
adjacencies flagged for an owner read, both no-collision: `2608.30238`
(unrestricted multiplicative complexity of 4-term polynomial multiplication,
MC = 9, Lean 4 — a different tensor and a different measure than
$M\langle3,3,3\rangle$) and `2608.30273` ($\Theta(C_7)\ge3.25883262$ — a
lower-bound-only explicit construction in kg's theta neighborhood with no
theta/SDP upper bound; $K_G$ untouched).

**Third coverage defect in our own sweep instrument, and the sharpest yet.** The
2026-08-30 sweep printed ONE `max published` (`2026-08-27T17:53:45Z`) for a
TWO-category sweep. Owner re-derivation: that is exactly `cs.CC`'s 08-27 stratum
maximum, while `cs.IT`'s 08-27 maximum is `2026-08-27T20:14:37Z`. Three `cs.IT`
submissions sat in that 2 h 21 min shadow — `2608.27565` (`18:00:04Z`, four
seconds past the 14:00 ET cutoff), `2608.27635` (`19:19:19Z`), `2608.27682`
(`20:14:37Z`) — and the Monday morning query then excluded the entire 08-27 day
through its `202608280000` lower bound. Verified property: all three are
simultaneously **after** the prior sweep's printed maximum and **below** the next
query's lower bound, so two instruments missed the same three items for two
different reasons. All three triage clean (surface-code neural BP decoder,
Siamese-net interference, bivariate bicycle codes — quantum/wireless, no Goppa,
no target), so no result is affected; the defect is coverage bookkeeping. Second
instrument finding: `eccc.weizmann.ac.il/report/2026/163` answers **HTTP 200**
with a `404 PAGE NOT FOUND` body, so a status-based absence probe would have
reported the report as existing; the content discriminator (`TR26-\d+`
identifiers present for 162, absent for 163) is what settles it. The agent's
report described the body correctly; the owner's first reading of it as a status
code was wrong and is corrected here. Both lessons are now rule 17d in
`README.md`. Sweep bookkeeping additions: ECCC publishes on some weekend days
(TR26-153, a Saturday), so weekday assumptions may never bound a sweep; and
weekend-DATED submissions exist while weekend ANNOUNCEMENTS do not — the sharper
phrasing replaces the 2026-08-30 entry's looser wording, which stands as written
per the running-tally rule.

Next window: Tuesday 2026-09-01, arXiv `cs.CC`/`cs.IT` bounded per category with
a printed boundary each, ECCC from ceiling TR26-162 bounded by publication date.

### Ledger-entry verification made mechanical before the next sweep, 2026-08-31

The gap was narrower than "unverified claims" and worse than it looked. Row
integrity was genuinely machine-checked all along — `--validate-resume`,
per-line hashes, the static guard, the checksum ledgers. But the **aggregate**
counts every ledger sentence cites were summed once in an ephemeral owner
session, written into `manifest_final.json` and `report.md`, and only then
hashed. Searching the runner-produced artifacts for those aggregates returns
zero hits in `summary.json`, `static_guard.json`, `run.log` and
`orbit_rows.jsonl`: the runner emits per-row components and a status/timing
summary, never the totals. So the checksum chain **froze the owner's
arithmetic without ever checking it**, and a derived file became an
independent authority for numbers whose only source is the canonical rows —
the SSOT inversion, not a trust failure.

[`tools/verify_campaign.py`](tools/verify_campaign.py) closes it, as owner-only
infrastructure that reads finished artifacts and is imported by no campaign
runner. It re-derives checksum ledgers, per-line row hashes and their key
bijection, declared source hashes, outward-decimal/binary-rational identity,
Arb-ball enclosure by the archived outward endpoints and width, CSV tables as
order-independent row-wise projections with every candidate source field named,
and the sums/maxima/histograms of every leaf shared by all rows. `--claims`
diffs stated numbers against derived ones and exits non-zero on disagreement.

Applied to the three current `delcap/` campaigns: `q3-invariance-correction`,
`q4-total-output-orbit-mass`, and `q3-n11-total-output-orbit-mass` all return
`PASS`. The 52 numbers this session published for `q=3,n=11` across
`RESULTS.md`, `PROGRESS.md` and `delcap/README.md`—row and status counts, the
four interval triples, 1,417,176 dual evaluations as 708,588 per candidate,
634,415,384 dual terms as 317,207,692 per candidate, 316,716,419 primal
conditional entries, 2,125,764 input and 9,565,920 output generator checks,
zero overlap failures on the primal and each dual separately, 828 skipped
zero-mass words, the bump counts 5,448/0/1/4, max RSS 929,677,312, the census
tuple, and three file hashes—were re-derived from the frozen artifacts and all
52 match
(`MACHINE-VERIFIED`). Claims file:
[`delcap/scratch/claims_q3_n11_ledger.json`](delcap/scratch/claims_q3_n11_ledger.json).

Negative controls, both on copies, never on frozen evidence: perturbed claims
(off-by-one count, wrong row total, wrong file hash) were each rejected with
the artifact value printed; one flipped byte inside a canonical row was caught
three times over, by the file ledger, the per-line hash, and the CSV
projection. The superseded 2026-08-30 `delcap/` campaign is correctly reported
as having no per-line ledger — file-level integrity only — which is a fact
about that older artifact, not a repair target.

Trimmed on review rather than defended: a `sums_by_leaf` namespace was cut. It
produced 74,720 leaf keys for a four-row campaign and merged distinct
quantities that share a leaf name — the single claim using it had silently
pooled the primal and both dual overlap counters, now three explicit paths. The
order-independent CSV matcher was kept because it is load-bearing on existing
evidence: `q3-invariance-correction`'s table is a genuine projection whose row
order differs from its canonical `.jsonl`, so positional comparison falsely
failed it. No `n=12` tooling was written; none is needed until that campaign
is frozen.

Standing protocol for the next multi-target sweep, unchanged from the launch
contract and now enforced: Main is the single writer of `README.md`,
`RESULTS.md`, `PROGRESS.md` and `docs/`; each agent writes only its own
`cs/<target>/`; sweep output lands non-authoritative; and each campaign's
numbers pass `verify_campaign.py --claims` before its ledger sentence is
written, one campaign at a time.

### `delcap/` q=3,n=11 extension completed; 4/4 strict finite-$n$ improvements, 2026-08-31

The named successor to the q=4 grid is closed. The frozen `08ecbc49…`
runner certified all four registered `q=3`, `n=11` rows at
`d in {1/2,1/5,1/10,1/20}`. Every row is `CERTIFIED`, every selected dual is
`ba_total_orbit_mass`, and every conservative Arb interval strictly beats
both published Tavakoli-Nguyen-Bose finite-$n$ endpoints. Outward width upper
bounds range from `1.7049647766059455e-06` at $d=1/20$ to
`7.068320598920295e-06` at $d=1/5$.

All four rows reproduce the exact census tuple
`(14884, 22450, 6148309, 3573542, 9721851)` and the independent Burnside
counts 14,884 input and 22,450 output orbits. The rows record 1,417,176
direct full-input dual evaluations, 634,415,384 positive-$W$/positive-$D$
dual term checks, 316,716,419 direct primal conditional entries, and
11,691,684 expanded input/output generator checks, with zero
representative/full-word overlap failures. Immutable-resume validation
passed with exactly four rows; owner checks re-derived every
newline-inclusive row hash, re-verified each outward decimal against its
archived exact binary rational, and confirmed the exact CSV projection
(`MACHINE-VERIFIED`).

Production took `1151.072 s` wall / `1150.401 s` CPU; maximum row-recorded
RSS was `929,677,312` bytes and no resource stop fired
(`COMPUTATIONAL-EVIDENCE`). Pre-statement and static evidence remain
byte-preserved; additive final manifest/checksum artifacts close the
campaign. Verification:
[`report.md`](delcap/campaigns/2026-08-31T20-34-52Z_04ed25bd-455f-4814-8288-5ee0106e0db8_q3-n11-total-output-orbit-mass/report.md).
Final checksum-ledger SHA-256 `0ba5629a…`; all 20 entries pass. Scope is
only the registered four-cell `q=3,n=11` box—no q>=4, `n!=11`, other-$d$, or
asymptotic claim. Named next campaign: `q=3,n=12`
total-output-orbit-mass extension under a separate prospective freeze.

### `delcap/` q=4 grid completed; 24/24 strict finite-$n$ improvements, 2026-08-31

The frozen total-output-orbit-mass runner completed all 24 registered cells
for `q=4`, `n=5..10`, and
`d in {1/2,1/5,1/10,1/20}`. Every row is `CERTIFIED`, every selected
dual is `ba_total_orbit_mass`, and every conservative Arb interval strictly
beats both published Tavakoli-Nguyen-Bose finite-$n$ endpoints. Outward width
upper bounds range from `6.165071504436472e-09` through
`2.3125135592425783e-05` bits/symbol.

The rows record 11,182,080 direct full-input dual evaluations,
3,585,643,216 positive-$W$/positive-$D$ dual term checks,
1,792,821,608 direct primal conditional entries, and 61,501,392 expanded
input/output generator checks, with zero representative/full-word overlap
failures. Immutable-resume validation passed with exactly 24 rows. Owner
checks verified every newline-inclusive row hash and the exact CSV projection
(`MACHINE-VERIFIED`).

Production took `4651.466 s` wall / `4649.094 s` CPU; maximum row-recorded
RSS was `1,441,169,408` bytes and no resource stop fired
(`COMPUTATIONAL-EVIDENCE`). Historical pre-production evidence remains
unchanged; additive final-manifest/checksum artifacts close the campaign.
Verification:
[`report.md`](delcap/campaigns/2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/report.md).
Final checksum-ledger SHA-256 `d2660432…`; all 23 entries pass. Scope is
only the registered 24-cell q=4 box—no q=3, q>=5, n>=11, other-$d$, or
asymptotic claim. Named next campaign: `q=3,n=11`
total-output-orbit-mass extension under a separate prospective freeze.

### `mm3/` new-decomposition landscape exhausted; minima 56 and 58, 2026-08-31

The frozen standard-action census reproduced all five startup totals
55/58/56/59/60 and every 729-entry Brent tensor before enumerating
`mws59` and `stapleton60`. It admitted 5,796 and 4,800 data triples,
corresponding to 1,922,973,696 and 1,592,524,800 valid orientations. Every
admitted triple received an exact integer lower bound and complete subset-DFS
floor decision.

`mws59` has certified range 56–74, upper median 67, two minimizing triples,
and nonmonomial-involving minimum 56. `stapleton60` has range 58–76, upper
median 67, nine minimizing triples, and nonmonomial-involving minimum 58.
There are zero lower bounds at most 54. Hence no new LB-55 triple enters from
either decomposition. `paper55` remains the only verified 55-addition circuit
in the expanded swept set, while `sun56` remains LB 55 / UB 56; this is not a
claim that `paper55` is the unique LB-55 data-triple minimizer.

An independent postcompute validator regenerated the pair masks, all survivor
sets, side bounds, totals, histograms, medians, minimizers, crosschecks, and
aggregate counts; 476 assertions pass. Its first attempt made a schema
assumption, aborted without a verdict, and left an empty write-once output;
the corrected retry and failed attempt are preserved. Precompute checksums
pass 17/17 and final checksums pass 31/31. Final-manifest SHA-256
`f33c4161…`. Scope is only the two named decompositions, ternary
determinant-`+/-1` `T^3`, all sigma powers, and the registered standard
action; no universal no-54 or upper-bound synthesis claim follows.

### `mceliece/` exact direct-route toy census exhausted; P2–P4 fail 24/24, 2026-08-31

The released scalar campaign passed every control and durably completed the
contiguous `beta=8,...,21` prefix before its registered 7,200-second wall cap
interrupted `beta=22`. Its 14-cell prefix is preserved as partial evidence;
the 22-row partial checksum ledger passes.

The separately preregistered exact table-indexed continuation then passed
26/26 controls and exhausted all 24 GF(32) cells in 1,622.94 seconds. Exact
RREF return tuples, nullspace basis order, and ranks matched the scalar
reference on 1,004 systems; a wrong multiplication table was detected; and
cached/uncached `beta=8` matched the parent semantic anchor. P1_U3 and P1_U4
pass in every cell and all 120 branches of each size. P2 fails 24/24 with 33
accepted labels against five expected; P3 fails 24/24 with 28 unmapped
labels; P4 fails 24/24 with `rank(E)=36`, `nullity=244`, and both admission
predicates false. Registered toy `genericity_pass` is therefore false.

Frozen C26 removed `complete` in addition to the two permitted provenance
fields. The deviation is preserved. A post-run corrective audit removes only
`header,t_utc` and compares all 14 parent-prefix cells, including
`complete=true`; every semantic payload matches exactly. This is not
retroactively called a pre-production control. Successor report SHA-256
`662e6024…`; state manifest SHA-256 `ac21cdf3…`; all static, state, final, and
parent-partial ledgers pass.
Exact census, predicate, control, and equality statements are
`MACHINE-VERIFIED`; wall times are `COMPUTATIONAL-EVIDENCE` only.
Verification:
[`report.json`](mceliece/campaigns/2026-08-31T16-12-46Z_vectorized-rref-replay/report.json).
Scope remains the one toy population outside
current ePrint 1786 Table 1 and Assumption 1; no paper-scale premise is
verified or refuted.

### `omega/`, `kg/`, and `delcap/` follow-on campaigns integrated, 2026-08-31

- **`omega/` two-rung replay:** all three corrected-core child commands
  completed, then the frozen aggregate launcher exited 1 on
  `KeyError: 'omega_cert_upper'`. The rung-1 schema instead has separate raw
  (`2.371340083602922`) and absorbed-slack (`2.3713411672715115`) values. A
  preregistered hash-gated recovery compared both against rung 2 with Lemma
  (`2.3715538358544617`): both orderings hold and all values are below
  `2.371866`. The failed key was not replaced; formal status remains
  **SUSPENDED** for the two original proof blockers. Recovered artifact
  SHA-256 `8160e846…`; all 22 checksum entries pass.
- **`kg/` midband refinement:** 6/12 tiles certified after 56,451 envelope
  evaluations. Five adjacent tiles extend the `[1.75,3.5]` certified prefix
  through `3.107728208020454953573248`; one terminal tile through `3.5` also
  certifies. Six tiles remain OPEN, leaving contiguous gap
  `[3.107728208020454953573248,3.49980671715929653996950394581]`. Result
  SHA-256 `f9760458…`; all nine checksum entries pass. Full band remains
  **PARTIAL**. Static source audit also retracts the proposed panel-wise `K_o`
  radius as a new mechanism: it was already present.
- **`delcap/` $q=4$ anchor:** the released `n=5,d=1/2` startup anchor passed
  with outward interval
  `[0.6664806108007938,0.6664806978587714]`, all registered group checks, and
  both dual candidates evaluated on all 1,024 inputs. Anchor SHA-256
  `c47aeb35…`; all 14 current checksum entries pass. It intentionally wrote no
  row artifact, so the full 24-row grid remains **COMPUTE PENDING**.

Process disclosures remain active: the high-band staging agent deleted three
agent-generated bytecode caches without required confirmation; they were not
restored and no claim depends on them. The earlier octonion agent's two
untracked-draft deletions were restored byte-exactly.

### `omega/` repaired candidate-ray replay complete; all five moves worsen, 2026-08-31

Campaign
`omega/campaigns/2026-08-31T14:29:40Z_OmegaGateCEntropySolver_v12_d2d3d4d5/`
replayed the base and all divisors `16,8,4,2,1` under the repaired outward
interval core. It hard-hashes local dependencies, converts exact finite Arb
dyadics without binary64 rounding, uses L-BFGS-B only as a warm start, and
accepts entropy multipliers only through an outward Arb KKT certificate.

All 1,152 entropy certificates passed; the largest KKT upper bound was
`1.4456822252392066e-10 < 1/2^32`. Every full-objective moved-minus-base
interval is strictly positive, from
`[6.475289389514407e-9,6.475289390521749e-9]` at `R/16` through
`[1.0365687869954472e-7,1.0365687870055208e-7]` at `R`. Thus the named ray
contains no descent at the five registered points. Status remains **OPEN**:
the ray preserves a mass strictly below one and one worsening ray proves
neither constrained local optimality nor an exponent bound. Result SHA-256:
`0403e652099a217d421c9bb1aa7fd56f07fb6daff72ddf14273d6a6dff063122`.

The v11 production is preserved as **ABORTED** at its stale old-objective raw
anchor. A fresh Fable audit of v12 reported no release blocker; owner checks
resolved its two residual dependency questions, and the 16-entry checksum
ledger passes. Only the candidate-ray statements are restored. The historical
two-rung and box enclosures remain **SUSPENDED** pending their own corrected
replay.

### `kg/` direct Lemma-D.4 restart partially re-certifies four bands, 2026-08-31

Campaign `kg/campaigns/20260831T082425Z_kg_direct_d4_restart/` removes the stale
startup anchor and applies the direct 256-bit Lemma-D.4 envelope. Bands
`[1,1.3]`, `[1.30,1.45]`, and `[1.45,1.75]` pass on all 14, 6, and 10 tiles,
with worst certified margins `+9.495739844715277e-7`,
`+1.95705024245827645e-6`, and `+3.10559377267758493e-6`.

Band `[1.75,3.5]` is **PARTIAL**: 25/36 tiles pass, covering a continuous
prefix through `2.81476518658164445646615405049` and the disjoint terminal
tile beginning `3.49980671715929653996950394581`. The 11 open interior tiles
hit depth 14 and 32,768 panels. Their negative upper-envelope lower margins
are failures to certify, not counterexamples; no refutation trigger fired.

The `[3.5,4.083]` baseline completed with 0/8 tiles certified after 3,900
envelope evaluations: all eight reached the panel cap. The closest open
upper-envelope lower margin is `-1.59862092695915505593006178343e-5`; the
most negative is `-1.12004131907617617264576323118e-4`. These remain
upper-envelope failures only, not refutations. Result SHA-256:
`35030d91f4a220f348991cd6731951118f8b7887d5bd47991e219b510980968e`.

The `[4.083,6]` panel-cache refinement completed with 0/20 tiles certified
after 14,047 evaluations. All reached depth 22 and the panel cap. The closest
open margin improved to `-1.97261882086860844646460017699e-7`, but no tile
closed; result SHA-256 is
`7e144344c11127ef3a069ca28053719b50ba1b9d730de616c22981c349f61b3c`.
This is again an upper-envelope failure, not a refutation.

The widened unresolved-interior refinement is now the sole low-priority,
thread-one CPU process. Its static refreeze corrected decimal `fmpq`
construction before release, passed its 12-tile cover and release-refusal
checks with zero cache entries/results, and passed its six-file checksum
ledger. Full Gate B remains open.

Process disclosure: the highband staging agent deleted three agent-generated
`__pycache__` files without the required confirmation. They were compiled
cache artifacts and were not restored; the deletion and all staging probes are
preserved in the campaign's additive correction. No claim depends on them.

### Held extension campaigns pass static gates only, 2026-08-31

The deletion-channel $q=4$ campaign
`delcap/campaigns/2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/`
passed a fresh no-iteration guard: the $q=4,n=5$ census anchor, 130,560 joint
equivariance checks, required breaking plants, and exclusion of the old
$q=3,n=11$ path. It created no anchor, BA iteration, or result row; production
remains held.

The McEliece exhaustive-genericity campaign
`mceliece/campaigns/2026-08-31T08-44Z_DA168C79/` was refrozen with explicit wall
caps and passed its released `--check-only` path plus all five static hashes.
No state directory, control result, C9 census, or production result exists.
The MM3 different-decomposition campaign
`mm3/campaigns/2026-08-31T083323Z_55447c11-5887-4fdd-aa2f-30e31f2332e8_b3046e6a6c50/`
passes all 17 precompute checksums and its released `--check-only` perimeter;
that mode imported no target module and created no result.

All three campaigns remain pending the exclusive CPU slot. These are
preflight facts, not new theorem or certificate claims.

### `omega/` outward-enclosure defect found; rigorous labels suspended, 2026-08-31

Independent frozen-source review found a load-bearing defect in the shared
endpoint-pair interval kernel: `from_two(a, b)` takes `b.lower()` as its upper
endpoint.  The owner confirmed the code path in the live and frozen copies.
Nonzero-radius Arb upper expressions can therefore be rounded inward.  The
two-rung endpoints, certified-box/full-box enclosures, and candidate-ray
differences are preserved as old-code outputs but **SUSPENDED as formal
interval certificates** pending corrected replay.  Float reproductions and
exact integer/rational checks are unaffected.

The v9 generalized-entropy run independently aborted at the base point on a
transposed feature matrix before its first certificate.  v10 repaired that
shape but was rejected without a production run: its binary64-only converter
would abort on the registered higher-precision dyadics, and its dependency
still carried the enclosure bug.  A clean v11 freeze is in progress with
`b.upper()`, exact Arb `man_exp()` extraction, true-Arb $1/e$, a production
`__debug__` gate, residual-system-aware cache keys, frozen code and data
inputs, and explicit regression plants.  No strict sign has been established.

### `oct-rank/` Route A closed at 12; Route F isolated to 13-or-14, 2026-08-31

The former commutator “16” route is now algebraically and source-theoretically
closed. Exact coefficient identities plus a human derivation prove

$$[L_{\bar u}L_v,L_{\bar u}L_w]^2
=-4N(u)\det\operatorname{Gram}(u,v,w)I,$$

so every independent real octonion triple has commutator rank eight. Direct
owner reading of Landsberg Theorem 6.1.1 confirms the factor-$1/2$ inequality:
rank eight implies only $r\ge12$. The naked $r\ge8+8=16$ paraphrase is false;
its quaternion analogue gives eight against known rank seven.

For the fixed $(1,i,j)$ triple, exact global conjugation gives two shared-first-
factor quaternion blocks. The known $\tau$ witness gives upper 14 and the
substitution/pencil chain gives lower 13. Applicable direct-sum or
multiplicativity machinery was not found; five CP and three
commuting-extension numerical searches missed but prove nothing. Correct
status: **OPEN in $\{13,14\}$**, with precisely a rank-13 witness or a rank-13
impossibility proof missing.

Campaign:
`oct-rank/campaigns/2026-08-31T08:02:18Z_routeAF/`; owner checksums accepted
23 payloads. The full octonion multiplication interval is unchanged at
$[18,25]$.

Two process breaches are part of the evidence record. The numerical searches
overlapped `kg/`'s exclusive slot for 105.30 s and 32.67 s, so no timing claim
is admissible. The agent also deleted two untracked drafts without the
required user confirmation; both were restored byte-exactly from the edit
transcript, marked `UNEXECUTED / NOT EVIDENCE`, and no lasting file loss was
found.

### `omega/` slope replay closed; direct candidate witness queued, 2026-08-31

The Gate-C diagnostic campaign repaired the endpoint replay and passed the
frozen V1 gate. All seven $R$ branches match, with
$R=2.8170035674609757$, $M=2.0942543887102634$, and raw
$\Omega=2.3715538358350807$; residual
$4.440892098500626\times10^{-16}$.

The preserved $+0.9008715$ failure is now retracted as nonreproducible. Its
advertised source crashes on undefined `s_const` before the first aggregate
checkpoint. The owner-suggested penalty explanation is also falsified:
zeroing the narrow Lemma-1 contribution moves the endpoint by only
$2.029648699330977\times10^{-6}$.

Before the decision pass, four independent bugs were removed: extra radius
scaling, wrong LP-dual dimensions, substitution of a dense SVD basis for the
exact integer kernel, and omission of the arbitrary-$y$ full-box instrument.
Binary64 inward endpoints and single-argmin branch pinning were replaced by
exact Arb endpoint formation and hulls of every possible minimizing branch.

The midpoint LP now supplies a numerical candidate
$-1.577332901516852\times10^{-7}$, but the rigorous full-box slope enclosure
is $[-5.841014555619999,5.841014555619999]\times10^{-6}$. It straddles zero:
**FAILURE TO CERTIFY / OPEN**, not local optimality and not an improving
witness. Frozen campaign:
`omega/campaigns/2026-08-31T07:57:00Z_OmegaGateCDiag_4859327/`; owner
checksums passed all 17 payloads.

The next bounded question is prepared as `OmegaGateCCandidateWitness`: lift
the single LP direction through the exact integer kernel, predeclare a dyadic
step ladder, and directly compare the full branch-safe Arb objective at the
center and moved points. It is held behind the single-heavy-process slot.

### `delcap/` symmetry-certificate defect found — 19 $q=3$ rows suspended, 2026-08-31

**A full-alphabet regression stopped the extension campaign before any new
$q=4$ row.** The old dual path rounded each entry of the output reference
distribution $D$ independently, then maximized KL divergence only on input
orbit representatives. At $(q,n,d)=(3,6,1/5)$, the exact denominator-$2^{30}$
snap split flat output orbit 51 (size 12): its entries contain both $578{,}028$
and $578{,}029$. The representative reduction is therefore inadmissible.

The failure is material, not cosmetic. A planted $(3,3,1/2)$ case gives a
strictly smaller representative-only maximum than the true full-alphabet dual.
The audit's first frozen row, $(3,6,1/2)$, did pass exact orbit equality,
$3{,}279$ generator checks, all $729$ input divergences, and equality of
representative/full maxima. The second row failed; the mandatory stop left
18 unrun. Status of the old 20-row ladder: **1 re-certified, 19 suspended**.

Evidence:
`delcap/campaigns/2026-08-31T08:01:53Z_9108c7d1-73f8-40dc-8c4a-4db1c6d9dcbf_q4-orbit/frozen_snap_invariance_audit_v3_failure.json`.
The target README and root results ledger now retract the unsupported
certification immediately; no Tavakoli–Nguyen–Bose claim is affected.

**Repair in progress:** snap total masses in exact **output-orbit**
coordinates and set $D_y=D(O)/|O|$, then prove positivity, normalization,
output-orbit invariance, and agreement between representative and every-input
dual maxima. Any symmetry-compressed primal gets an independent
full-alphabet equality check. All prior 20 rows must pass this corrected path
before the held $q=4$ extension can start.

### `delcap/` corrected replay complete — 20/20 $q=3$ rows restored, 2026-08-31

The repair above has landed as a replacement campaign:
`delcap/campaigns/2026-08-31T09:13:25Z_7dc5babe-5e1b-44fa-a9e5-03d2f10183b6_q3-invariance-correction/`.
It snaps total input- and output-orbit masses exactly, proves normalization
and $G$-invariance, and compares the compressed primal and dual against direct
full-alphabet evaluations. All 20 fixed rows at
$q=3$, $n=6,\ldots,10$, and
$d\in\{1/2,1/5,1/10,1/20\}$ now certify strict improvements over both
published sandwich endpoints.

The replay found 16/20 old `ba_word` candidates non-invariant; the corrected
candidate is `ba_total_orbit_mass` in every row. It checked 1,058,508 direct
input KL values, 5,821,704 exact group actions, and 237,744,816
positive-$W$/positive-$D$ terms, skipping no positive-mass input. Owner
`sha256sum -c` verification passed all 17 payloads. The certified outward
widths span $1.2226475645090283\times10^{-8}$ through
$2.220013228610702\times10^{-5}$.

Status: the old 19-row suspension is closed by new certificates; the
defective endpoints remain retracted. The $q=4$ extension is separately
frozen but **NOT RUN**, and no asymptotic or published-theorem claim changes.

### `rs-pe3d/` H-GATE mechanism proved — exact distance of the lifted-line sum code, 2026-08-31

**Correct theorem (HUMAN-AUDITED, machine-anchored):** let
$A_i=\mathbb F_q^{s_i}$, let $C_i\subseteq A_i$ be the stated
$\Lambda_i$-diagonal generalized Reed–Solomon code of dimension $t_i$, and let

$$V=\sum_i A_0\otimes\cdots\otimes C_i\otimes\cdots\otimes A_2.$$

For distinct evaluation points, nonzero diagonal entries, and
$1\le t_i\le s_i$,

$$d_{\rm Ham}(V)=\min_i d(C_i)=\min_i(s_i-t_i+1).$$

Equivalently, **every** point support $S$ with $|S|<d_{\rm Ham}(V)$ has
$V\cap\mathbb F_q^S=\{0\}$, while **some** line support of size
$d_{\rm Ham}(V)$ contains a nonzero word. Therefore the union window over all
supports of weight at most $w$ is empty iff $w<d_{\rm Ham}(V)$. This is not an
individual-support biconditional for $|S|\ge d$.

**Proof mechanism, owner-derived and independently re-derived by the target
agent.** With $Q_i=A_i/C_i$, right exactness gives
$V=\ker(\pi_0\otimes\pi_1\otimes\pi_2)$. Any fewer than $d(C_i)$ coordinate
cosets in $Q_i$ are independent, since an arbitrary relation would produce a
nonzero $C_i$ word on smaller support. For a target point in a support of size
$w<\min_i d(C_i)$, assign every other point to one coordinate where it differs;
on the resulting **sets** of coordinate cosets, quotient-dual functionals
separate the target. Their tensor pullback annihilates $V$ and isolates the
target coefficient, forcing the supported word to vanish. Sharpness uses a
degree-$(t_i-1)$ GRS polynomial vanishing at exactly $t_i-1$ evaluation
points, lifted on one line.

**H-GATE is now resolved on its actual slice.** At $t=(1,1,1)$,
$d_{\rm Ham}(V)=\min_i s_i$: the weight-$\le3$ window is empty exactly when
$\min_i s_i\ge4$, and a line word makes it nonempty when $\min_i s_i\le3$.
For general positive $t$, the predicate is
$\min_i(s_i-t_i+1)$, **not** $\min(s)$. Pairwise coprimality and nonidentity
$\Lambda$ do not affect this distance statement. The $t_i=0$ zero-code case
is separate.

**Exact controls reported by the frozen campaign:** dimensions $13,22,37$;
quotient-kernel dimension $37=64-3^3$; complete below-distance emptiness over
$16$, $30$, and $1{,}830$ supports; $82{,}160$ weight-3 supports empty at
$(41,(4,4,5))$; $5{,}050$ weight-$\le2$ supports empty at
$(61,(4,5,5))$; all $14/14$ and $41/41$ corner supports nonempty when an
order is one; prior nonempty anchors $3486=3486$ and $917=917$; exact line-word
$\delta$ controls $4$ and $2$. Frozen at
`rs-pe3d/campaigns/2026-08-31T08-26-54Z_hgateMechanism/`; owner
`sha256sum -c` accepted all four checksummed files. The consolidated replay is
CPU-bound and remains queued behind the repository's single-heavy-process
slot.

**The frozen prose must be cited only through the target README's owner
correction.** Owner review found seven errors: per-support “iff” instead of
union-window semantics; `1_T` instead of an arbitrary relation vector;
$s_i-1$ instead of $t_i-1$ roots; $N-\prod t_i$ instead of
$N-\prod(s_i-t_i)$; a false uniqueness claim for the triple-sum
representation; no separate git commit proving the claimed preregistration;
and PP2 mistyped as $(2,3,4)$ instead of the actually replayed $(2,3,7)$.
The immutable campaign preserves every original; the appended README
corrections supersede them. One 52.9-second census also overlapped `kg/`'s
heavy slot, a recorded resource-contract breach; its wall time is not
performance evidence.

**Rule-7 boundary:** this proves the support-distance mechanism, not
Conjecture 4.2, not any global $\rho_{\rm inst}$, and not an optimal ratio in
the weight-$\ge4$ window. At $\min(s)\ge4$ and $t=(1,1,1)$,
$\rho^{\rm window}$ for the old weight-$\le3$ program is **undefined**, not
zero. The first nonempty window starts at weight $d$; optimizing there remains
OPEN.

### `mceliece/` live-source event — ePrint 1786 changed route four revisions after cache, 2026-08-31

**The source audit changed the question before the new campaign began.** At 07:56 UTC, live
ePrint `2026/1786` was **340,507 bytes**, sha256
`12e244c760a068d74bc2784c79cfe0b0a6eef698d7bf10994fadb5c7bc50731d`; cached
`mceliece/prim/pe1786.pdf` was **556,874 bytes**, sha256
`ee0c447b251d5c18735f4687ac8b476e9aeb5dcd57c5660a4fc3fb57e4da9d21`.
The archive identifies the cached bytes exactly as version `20260828:000118`. Four PDF
updates then landed on 30 August, culminating at `20260830:172206`. This is not packaging:
the current paper's own changelog says Sections 2–5 replaced the singleton tangent/anchor route
with a **one-chart direct locator-recovery route**, removed the superseded singleton tables,
and expanded Assumption 1.

**Current headline, source-conditional:** complete arithmetic/state estimates
$2^{126.77}$, $2^{145.22}$, $2^{137.48}$, $2^{136.65}$, $2^{137.48}$ bit operations
for the five Classic McEliece parameter sets, each below its cited NIST classical-gate
reference. The paper explicitly conditions them on four finite-instance premises:
canonical rank-two branch; exactly the common Frobenius branches in the preliminary test;
one globally compatible projective label per coordinate after four holdouts; and full plus
block-projected kernel-rank admission. Failure of the third/fourth on $r$ unavoidable
coordinates produces the paper's explicit fallback obstruction $N_{\rm direct}\ge2^r$.
**No target-scale public relations, locators, Goppa polynomial, or key are reported.**

**Arithmetic replay, owner-independent.** The current `mjosaarinen/mccost` artifact's exact
dimensions, integer costs, and five printed totals replayed. Its threshold path uses binary64,
so the owner separately rebuilt the aggregate solver-failure expression as an exact rational
and evaluated only its display logarithm in 300-bit outward Arb. Current solver degrees all
strictly pass $2^{-128}$, while the immediately smaller admissible extension degree strictly
fails:

| set | degree | certified $\log_2 p_{\rm fail}$ | predecessor | certified predecessor |
|---|---:|---:|---:|---:|
| 348864 | 240 | $-135.1793047803\ldots$ | 228 | $-123.1793047803\ldots$ |
| 460896 | 260 | $-138.1954385517\ldots$ | 247 | $-125.1954385517\ldots$ |
| 6688128 | 247 | $-132.5182649166\ldots$ | 234 | $-119.5182649166\ldots$ |
| 6960119 | 247 | $-133.2543840595\ldots$ | 234 | $-120.2543840595\ldots$ |
| 8192128 | 247 | $-132.5182649166\ldots$ | 234 | $-119.5182649166\ldots$ |

Evidence split: solver-threshold arithmetic **MACHINE-VERIFIED**; the map from the paper's
structural premises to those dimensions and costs remains **CITED-DEPENDENCY, CONDITIONAL**.
Reproducing the cost program is not evidence for Assumption 1.

**The old gate is semantically disjoint.** Existing gate B measures Apon Lemma 9's
$\Delta_{p,q}\ne0$ Wronskian condition. It tests none of current 1786 Assumption
1(1)–(4). Same subject area is not semantic overlap. The new campaign was stopped before
computation and redirected to provenance repair plus the current four-holdout/rank premises.

**Second provenance defect:** `pe1747.pdf` on disk is the old 245,667-byte copy
(`e9fc4fff…`), while the file named `pe1747_v3_2026-08-30.pdf` is the current
264,802-byte copy (`f3e83fe7…`). The README describes the opposite canonical/superseded
roles. Both byte versions are retained; the target owner was ordered to repair naming
non-destructively and never delete either content version.

**Freshness controls:** a single arXiv `id_list` query returned the current versions of 14
target and adjacent primary papers; none has an `updated` timestamp after its recorded read.
Remote/current ePrint `1810`, `1778`, and `1630` are byte-identical to their cached PDFs.
ECCC `TR26-150` is still 1,119,003 bytes with sha256
`fb51d7dc84860f115d88ec47e92799e9c949ac90cbe714e8df7d3904c930af87`,
matching the target README's recorded prefix. These negative controls matter: the audit detects
both equality and change, rather than declaring every fetched source “fresh.”

### Frontier forward sweep — Friday delta announced Monday, zero target collisions, 2026-08-31

**Coverage coordinates named exactly.** At 07:51 UTC, arXiv API queries over
`submittedDate:[202608280000 TO 202608312359]` returned **3** `cs.CC` entries and **14** `cs.IT`
entries, with one cross-list, hence **16 unique submissions** after the prior 2026-08-27 ceiling.
The maximum `published` timestamp was printed in both feeds: `cs.CC`
**2026-08-28T17:14:45Z**, `cs.IT` **2026-08-28T17:35:09Z**. There are zero submissions dated
29–31 August as of the sweep. These are Friday submissions visible in Monday's window, not
Monday-dated submissions. A separate `lastUpdatedDate` query returned the same all-v1 set and
no revision-only record.

**ECCC bounded by publication date, not paper date or report number.** The latest-reports feed's
`lastBuildDate` and newest `pubDate` are both **2026-08-30 18:05:01 +0300**, still `TR26-162`;
the full 2026 year index independently tops at `TR26-162`. Therefore **zero new ECCC publications**
after the existing ceiling as of this sweep. This wording matters: ECCC report numbers are not
chronological, and `TR26-162` itself has a paper date one day earlier than its publication date.

**All 16 arXiv records triaged against all seven targets — no collision, no gate affected, no
target opened.** Exact delta:
`2608.28556`, `2608.28539`, `2608.28500`, `2608.28469`, `2608.28441`, `2608.28354`,
`2608.28222`, `2608.28179`, `2608.28133`, `2608.28086`, `2608.28085`, `2608.28068`,
`2608.27936`, `2608.27909`, `2608.27859`, `2608.27829`.

The closest title-level adjacencies were checked rather than discarded by keyword:
- `2608.28556`, tensorizable $f$-divergences, and `2608.28068`, finite-sample composite
  hypothesis testing, concern general divergence/testing structure — **not** deletion-channel
  capacity or any finite deletion channel, so `delcap/` is untouched.
- `2608.28222` concerns self-orthogonal/LCD embeddings over
  $\mathbb F_q+u\mathbb F_q$ — **not** binary Goppa full-support degeneracy, so
  `mceliece/` is untouched.
- `2608.27859`, 3-restricted matching-vector families, is the already-audited **ECCC
  `TR26-156` under a second identifier**. Matching the objects by exact title and author prevents
  a double-count; its LDC/MVF route does not touch `rs-pe3d/` product expansion.
- `2608.28354` is geometric discrepancy, not Grothendieck discrepancy; `2608.28133` is a
  Lean-checked theorem on bent partitions, not an open certificate gate; none of the remaining
  applied communications, quantum commitment, random-lift, or entropy papers touches
  `kg/`, `omega/`, `mm3/`, `oct-rank/`, or the other live gates.

**Three deferred “best-fit” candidates were then full-read and all rejected under the
repository charter.** `TR26-157`'s “tight” SCS inequality is already closed by its Theorems
3.1/3.3, including an explicit family proving every smaller endpoint impossible; its only
adjacent open problem is the 38-year-old greedy ratio-2 conjecture, with no useful failure
certificate. `TR26-159` publishes asymptotic sorting bounds and asymptotic open questions,
not a finite constant or table a Bellman search could improve. `TR26-156` has a certifiable
small matching-vector optimization, but feasible $m=6$ cells stop near $n=4$–5 while its
theorem starts at $n\ge36$ (and its predecessor at $n>2160$); no named finite problem,
published finite record, or load-bearing computation exists. The apparent LDC link was also
too loose: Efremenko's record-driving 3-query example is 4-restricted and uses a sparse
decoder. Exact corrections are retained inline in the scan rather than replacing the
original premises. **No decorative eighth target opens.**

The ceiling now is **submission date 2026-08-28 17:35:09Z for the union of `cs.CC`/`cs.IT`**
and **ECCC publication date 2026-08-30 18:05:01 +0300**. They are deliberately different
coordinates and are not normalized into a fictitious common “latest date.”

### `rs-pe3d/` P2 battery COMPLETE — one owner hypothesis killed, one confirmed prospectively, 2026-08-30

19/19 pre-registered instances adjudicated, `pre_statement.md` committed as the campaign's first file.

**Pattern P survives, not established.** PP arm (strong, predicts *exactly* $3/5$): **9/9**, primes to
601, $k$ to 32, $\beta$ to 16, every one with $\mathrm{wt}=3,\delta=5$ and per-witness exhaustive
closure. PN arm (weak, predicts merely "not $3/5$"): 7/7 landed, with PN8/PN9 giving $1/2$. Total
out-of-sample: **B11 + 16 rows, 0 failures.** The PN direction remains thin, which for a biconditional
is where it usually dies. The $(1/5)^9$ figure is labelled illustrative and explicitly not a p-value.

**H-COP FALSIFIED — the owner's own coprimality observation, reframed and killed. Record 0-for-13.**
H-COP.a holds; H-COP.b fails on $\{3/5,1,1/2\}$ — two distinct sub-1 values against a predicted
$\ge3$. A near-miss that would have been easy to write up as "broadly consistent"; the agent didn't,
and scoped the falsification to the landed set rather than overclaiming.

**H-GATE PROSPECTIVELY CONFIRMED — the first owner hypothesis to survive a prospective test today,
and it counts only because of the agent's procedure.** The weight-$\le3$ window is non-empty **iff
$\min(s)\le3$**, zero exceptions across every row ever run (owner-verified: non-empty $\min(s)\in\{2,3\}$,
empty $\min(s)\in\{4,5\}$). **PN10 was a committed prediction — addendum and order proof written at
49% of the census with zero hits so far** — and returned EMPTY as predicted. Written after the census
it would have been worth nothing. Owner verified all three emptiness censuses reproduce exactly as
$\binom{L}{1}+\binom{L}{2}+\binom{L}{3}$: $457{,}450$ ($L{=}140$), $288{,}100$ ($L{=}120$),
$3{,}658{,}900$ ($L{=}280$). The magnitude confound is **broken by the data**: $(3,5,16)$, product 240,
non-empty; $(4,5,7)$, product 140, empty.

**Agent-supplied precision that narrows the claim, recorded because it weakens it.** The crisp iff-form crystallized only *after* PN4/PN5/PN6/PN7 had landed, so **H-GATE's prospective content is exactly and only the PN10 read** — one committed data point, not twenty. Those four rows are labelled **retrospective** in the target README with the boundary kept explicit. The 20-row consistency is real but retrofitted. **An agent volunteering the distinction between its prospective and retrofitted evidence, unprompted and against its own headline, is the single behaviour this repository most depends on.**

**Recorded at exactly the strength earned:** one committed prospective test — the same status pattern P
had after B11 — **not** established; mechanism from the dual product codes OPEN; $s_i=1$ untested.

**The most durable output was not the assignment.** H-GATE tells this target what its reported domain
**is**: $\{$instances with an order $\le3\}$. Every ratio ever reported here, **including the frozen
13**, is conditional on that predicate, and nobody had stated it before today. The agent's own framing
names the next campaign precisely: weight $\ge4$ supports, *"where the $\min(s)\ge4$ instances actually
live."*

**Owner checks OC1/OC2/OC4 all confirmed** from the frozen witnesses rather than payload strings, with
the clarification that B04-B06 are **pairwise** non-coprime. **Two of the agent's own stale
post-freeze reporting strings disclosed and corrected as rule-5 addenda**, frozen originals untouched;
and the vectorized empty-window kernel was validated exactly against the frozen census code on
non-empty instances ($3486=3486$, $917=917$) **before** being depended on.

### `omega/` gate C — second FAILURE TO CERTIFY, three defects fixed, and an owner arithmetic catch, 2026-08-30

**Outcome: V1 FAILED again, gate C stays PARTIAL, the 21-dimensional question is UNCHANGED and
OPEN, and the index gains no certified number** — the honest consequence. The agent stopped before
the LP per its pre-registration rather than proceeding on an unvalidated aggregation, the second
time this target has honoured that gate.

**Real movement, though.** The predecessor's named repair was correct and is now closed: `M_low`
restored $0.14985\to2.09425439$ against the frozen $2.09425439$ at $8\times10^{-17}$, with level-2
parity to `vxxz24_float.PartLv2.evaluate_post` at $6.9\times10^{-18}$ across 1377 parts and a working
defect-detection control. **Two further predecessor defects were found and fixed en route, neither in
the handoff note** — a part-penalty `lam_sum` sign flip and an endgame endpoint-slope Leibniz error.
**Three defects in one aggregation path**, which says something about the path rather than about any
one agent. V1 residual fell from $+39.6246462$ to $+0.9008715$, a **43.98$\times$** reduction.

**Owner catch: the named next blocker cannot explain the residual.** The agent names the glob-region
Lemma-1 $\epsilon$ aggregation — frozen tree $\approx1.8\times10^{-15}$, replay $\approx3.08\times10^{-6}$,
a $1.71\times10^9$ discrepancy. That is certainly a defect, but its **absolute** size is
$3.08\times10^{-6}$ against a residual of $0.9008715$: it would need **$2.9\times10^5$** downstream
amplification to account for it. So either $\epsilon$ is genuinely amplified that far through a
Lemma-1 denominator or logarithm, or **a fourth defect remains unidentified**, or both. Fixing
$\epsilon$ under the second reading moves V1 by $\approx3\times10^{-6}$ and ends in an identical third
FAILURE TO CERTIFY. **The agent was instructed to settle which before writing any more code.**

**The control that makes the diagnosis airtight, promoted from the agent's evidence list.** It
reproduced the **frozen protocol standalone** at $2.3715538358351957$ against frozen raw
$2.3715538358350803$ — delta $1.155\times10^{-13}$, **inside** the $5\times10^{-13}$ V1 threshold. So
the frozen protocol **passes the very gate the replay fails**, on the same machine, data and anchor.
That rules out the environment, the frozen data, and the anchor value, and isolates the defect to the
slope replay path.

**Owner refinement, derived from the agent's own parity data rather than intuition — and it
reframes the next cycle.** The agent reports L2 R branch, `num_blocks` and `mat_size` all
**parity-exact**, with divergence isolated to the "Lemma-1 $\epsilon$/**penalty** aggregation path."
That path has **two** components and both the agent and the owner conflated them when naming the
blocker: $\epsilon$ is $O(3\times10^{-6})$ and needs $2.9\times10^5$ amplification to carry the
residual, whereas the `lam_sum` **penalty** is naturally $O(1)$ — and the agent had already found and
fixed **one sign flip in exactly that term**. A penalty mis-aggregation of $0.9009$ requires **no
amplification at all**. So the likely answer is **(c) both, with the penalty carrying essentially all
of the residual and $\epsilon$ a separate smaller bug.** One-pass test handed over: zero the penalty
contribution and re-run V1 — if the residual collapses toward $3\times10^{-6}$ the penalty is the
blocker; if it barely moves, the refinement is wrong and the stage bisection is needed.

**Also stated for the successor:** the residual is now $38\%$ of the raw value, not a rounding-scale
amount. A $38\%$ error is a **structural term missing or double-counted**, consistent with the three
defects already found and inconsistent with a last remaining $3\times10^{-6}$ term. Recommended
diagnostic: bisect over **stages**, substituting the frozen tree's own intermediates one stage at a
time until the residual collapses — a dozen stages rather than a search over 1377 parts.

### `omega/` gate C relaunched on a one-term repair, 2026-08-30

`OmegaSlope` left the cleanest possible handoff: machinery built and running, one named missing
term, and an anchor tolerance. Relaunched (`OmegaGateCFinish`) on exactly that scope — build the
level-2 `PartLv2 mat_size_contribution` in slope form (target `M_low` $\approx2.09425$, currently
collapsed to $0.14985$), pass **V1** against the frozen **raw** endpoint $2.3715538358350803$ within
$5\times10^{-13}$, and only then run the LP/remainder stages **already coded** in
`gate_c_endgame.py` to decide the 21-dimensional nullspace question.

**One instruction added from today's cross-target pattern:** derive the block from
`vxxz24_float.PartLv2.evaluate_post` **first-hand**, not from the handoff note's paraphrase. Two
defects today were caused by implementing from a predecessor's paraphrase — `kg` lost the paper's
2-D even-box structure that way, and `mm3` inherited a map applied outside its documented domain.
The handoff note is a pointer to the primary source, not a substitute for it.

**Pre-stated acceptable outcome:** a second honest FAILURE TO CERTIFY is fine; a certified-looking
claim resting on an unvalidated aggregation is not. And even total success caps **12.8$\times$**
short of the published $2.37155181$, so no record claim can arise from this route in any branch.

### `mceliece/` gate B FROZEN — the first measured bounds on Apon's own stated hole, 2026-08-30

**Zero degeneracy events across ~570 exact $\Delta$-verdicts.** Four claims at four distinct
strengths, kept separate as instructed:

1. **EXHAUSTIVE THEOREM at $(m,t)=(6,2)$:** 0/2016 over **every** monic irreducible degree-2 Goppa
   polynomial over $\mathbb F_{64}$ at full support. Owner verified the population is exactly
   $(q^2-q)/2=2016$, so this is the complete cell — **`MACHINE-VERIFIED`, no confidence interval**,
   since a CI on an exhausted population is a category error.
2. **Sampled layer:** 0/250 (10 cells $\times$ seeds 1000-1024), exact one-sided 95% Clopper-Pearson
   $\le0.0119$; owner reproduced $1-0.05^{1/250}=0.011911$, plus $0.013053$ at 228 and $0.023312$
   at 127.
3. **Adversarial special-$G$ families:** 0/127 ($42+21+64$, owner-checked), CP95 $\le0.0233$.
4. **NO forcing configuration found**, so **Apon's §3.6 hole is BOUNDED, not closed.** Stated as
   plainly as the positives. Neither disputant in the four-paper dispute had any of these.

**What makes the zero credible rather than an absence of evidence.** The agent built a **new**
instrument (rank-scan + completeness lemma $\deg_x\Delta\le D-1<n$) instead of reusing the inherited
route, and tested it against **5 counterfactual PLANTS** plus 3000 property cases against
brute-Wronskian ground truth — **catching a real iff-bug in its own first cascade draft before any
campaign number existed.** A predicate that always answered "nondegenerate" would have produced
exactly this headline; planting falsehoods is the only thing that separates the two. It also rebuilt
all 13 frozen instances **byte-exactly from seeds**, and correctly ruled $\Psi=C\!\cdot\!J$ off the
verdict path.

**It logged an incident against itself that nobody would have found:** an accidental touch on a
frozen campaign file, reverted and sha256-verified byte-identical.

**Two owner sharpenings, both of the agent's own results, both accepted and appended (never edited
through — rule 5 plus append-only).** (i) The $k=1$ cell was reported unreachable *"for $m\le11$"*;
in fact **no $m>1$ divides $2^m-1$**, so it is unreachable **unconditionally**. Canonical proof now
in the target README: take $p$ the **smallest** prime dividing $m$; $\mathrm{ord}_p(2)$ divides both
$m$ and $p-1$, but $\mathrm{ord}_p(2)\le p-1<p$ while every prime factor of $m$ is $\ge p$, forcing
$\mathrm{ord}_p(2)=1$ and $p\mid1$. Owner also verified numerically to $m=20{,}000$. (ii) The agent's
**own invariance theorem** — verdict is a pure function of $(m,t,G)$ — sharpens its **own** sampled
layer: the seeds vary neither support, basis, nor column order, only $G$, so the CP bound bounds the
rate over the **within-cell distribution of irreducible Goppa polynomials**, not over some unnamed
"sampling distribution."

**And the agent improved on the owner's question rather than answering it.** Asked whether the $t=1$
argument could be sharpened the same way, it distinguished: **unconditional at full support** (the
root of $Z+\beta$ is $\beta$, always in $E$), merely **rare at sub-support**, which is outside its
swept scope — so no sharpening is available there without a new sweep. A precise scope statement in
place of an overreach.

### In-flight status of the last two campaigns, 2026-08-30 — with two owner interventions

**`mceliece/` gate B (Apon $\Delta_{p,q}\ne0$ genericity rate).** ~570 exact $\Delta$-verdicts,
**zero degeneracy events anywhere.** Three layers kept separate as instructed: sampled 228/228
nondegenerate over 24 $(m,t)$ cells; **adversarial E1 EXHAUSTIVE — all 2016 irreducible $G$ at
$(m,q,t)=(6,64,2)$, 0/2016**; special-$G$ families 0/127. No forcing configuration found, so
**Apon's §3.6 hole remains OPEN at the scale measured** — bounded, not closed.

**Owner intervention 1 — E1 is a theorem, not a rate.** The count of monic irreducible degree-2
polynomials over $\mathbb F_{64}$ is $(q^2-q)/2=(4096-64)/2=\mathbf{2016}$ exactly, which the owner
verified. So E1 is the **complete population** of that cell, not a sample: for $(m,t)=(6,2)$ the
genericity condition holds for **every** irreducible Goppa polynomial. That is a finite exhaustive
theorem and **a confidence interval on it would be a category error.** Owner supplied exact
one-sided 95% Clopper-Pearson bounds for the two genuinely sampled layers instead: $N=228\to0.01305$,
$N=127\to0.02331$, with the instruction that these bound the rate over the **sampling distribution**,
not over the space.

**The agent's inherited-quantity disclosure was the best received today**, and one item is a finding
about its own work: it built a **new** verdict instrument (rank-scan + completeness lemma
$\deg_x\Delta\le D-1<n$), self-tested on 3000 property cases **and 5 counterfactual PLANTS** against
brute-Wronskian ground truth — and **it caught a real iff-bug in its own first cascade draft before
any campaign number existed.** Counterfactual plants are the control most skip: a test that only
sees true instances cannot distinguish a correct predicate from one that always returns True. It
also rebuilt all 13 frozen instances byte-exactly from seeds, proved *and* machine-checked
ordering/basis invariance, and correctly identified that $\Psi=C\!\cdot\!J$ belongs to the waterfall
census and is **not** on the gate-B verdict path.

**`rs-pe3d/` pattern-P battery.** `pre_statement.md` committed as the campaign's **first** file with
all 19 instances and their predictions fixed in advance, plus the owner's coprimality observation
pre-registered separately as H-COP. 12 of 19 landed exact: **PP arm 9/9** (all gave exactly $3/5$,
new primes $q=43,67,97,109,157,193,601$ up to $\beta=16$), **PN arm 3/3 landed**, zero misses.

**Owner intervention 2 — the tally must not be reported as "12-for-12", and this cuts against the
result.** The two arms are not equally informative. The **PP arm predicts one specific value** out of
the five observed across the corpus, so 9/9 there is the substantive evidence (crude uniform null
$(1/5)^9\approx5\times10^{-7}$). The **PN arm predicts "anything except $3/5$"** — a nearly free hit —
**and it is the under-tested arm**, 3 of 10 landed. For a biconditional, the weak direction is where
it usually dies. Required framing: *"9/9 on the discriminating arm, 3/10 on the weak arm with the
remainder outside the window."*

**And the battery's incidental finding was promoted above its own headline.** PN4 $(421,(4,5,7))$
has an **empty** weight-$\le3$ window — owner verified the exhaustive census size exactly:
$457{,}450=\binom{140}{1}+\binom{140}{2}+\binom{140}{3}$, so 140 lines, with sampled minimum weight
$\approx134$. If larger $s$ systematically pushes the dual product code's minimum weight above 3,
**the window cannot see those instances at all** — a rule-15 named-domain fact about the **entire**
gate-B programme, not just this battery: every ratio this target has reported, **including the frozen
13 rows**, is conditional on a non-empty low-weight window, and that condition correlates with small
$s$. More durable than another confirming row.

### `mm3/` off-diagonal programme COMPLETE — the "wall" was never there, 2026-08-30

**Result: 3,387,432,960 valid off-diagonal orientations over the named set, ALL certified $\ge55$,
ZERO at $\le54$, and the published 55-addition orientation is the UNIQUE minimizer of the entire
landscape.** Owner re-verified every figure.

| | `paper55` | `sun56` |
|---|---|---|
| survivor data-triples / $145^3$ | 5,700 (0.187%) | 4,510 (0.148%) |
| minimum certified total | **55** | **55** |
| triples attaining 55 | exactly **1** (all-monomial) | exactly **1** (all-monomial) |
| minimum among the rest | **56** | **56** |
| median | 66 | 68 |
| decided in | 929 s | 809 s |

Histograms sum exactly (5,700 / 4,510). Named-set total 3,387,432,960 instances — a **5,105$\times$**
extension of the inherited monomial-only surface $663{,}552$ (owner-verified:
$663{,}552\times5105$ reproduces it exactly). **All 10,208 non-monomial-involving data-triples
certify $\ge56$.**

**Four owner corrections, all from this one agent, all right.** (i) *"transpose used where inverse
was required"* was the symptom; the cause is **a formula that is a tensor automorphism exactly on the
orthogonal subgroup, applied outside it.** (ii) The data-node count is **145** ($=6960/48$, full
monomial group $B_3$), not the owner's 1160 (permutations only). (iii) The unswept off-diagonal
figure was **per-decomposition**, understating the region by exactly $2\times$. (iv) The owner's
per-factor collapse hypothesis was wrong by a factor of **exactly 5,700** — the survivor count it
had failed to see, since the prediction was precisely the single all-monomial triple.

**Two self-caught failures that carried the campaign.** The **A2 control caught the agent's own first
derivation** — the naive cyclic conjugation, which assumed the trilinear form is $\mathrm{tr}(UVW)$
rather than $\langle UV,W\rangle$ — failing Brent 10/10; compute halted per pre-statement and nothing
downstream touched the wrong map. Then a **bookkeeping false alarm** fired the $\le54$ trigger on
side-sums missing the $+14$ gap, and **only the anchor could distinguish it from a discovery**: a
stream of sub-54 totals is exactly what a real find looks like, and the all-monomial triple printing
41 instead of 55 falsified the field. Retracted inline before any external report.

**A full day was spent on a phantom disagreement** — owner and agent had the same map and the same
arrays throughout, differing only in whether the letters were written in standard or transposed-$L$
convention. Both independently enumerated all $4^6=4096$ assignments and found the **same unique**
all-ternary + Brent-passing family. **The lesson recorded with the formula in the index: state the
convention alongside any map, or two correct implementations will read as a contradiction.**

**Also overturned:** the record attack's *"the 48 survivors are exactly the monomials"* is
map-relative and false for `sun56` under the true family (576 valid, 528 non-monomial); its
**conclusion** (no $\le54$) survives and is now vastly extended, while its **stated reason** was
wrong. The 288 monomial decisions stood throughout.

### Three campaigns landed 2026-08-30 — delcap resolved, mm3 off-diagonal counted, oct-rank S3 filed

**`delcap/` Open Item 1 RESOLVED — 20/20 rows improve BOTH ends of the published TNB sandwich.**
Owner verified from `orbit_rows.jsonl`, not the report: `cert_lo` $>\max(\texttt{lb1},\texttt{lbplus})$
and `cert_hi` $<\texttt{ub}$ on **20/20 each side**, matching the artifact's own verdict field.
Widths $9.862\times10^{-8}$ to $1.320\times10^{-4}$ per symbol, all positive. Method: sparse
$G=S_q\times C_2$ orbit certificates on **symbol values** — the live reduction, separated from the
falsified position/type one by owner-run measurement in both directions. Three bugs fixed inline;
the important one produced a **NEGATIVE-width certificate** (dual below primal, sub-stochastic
per-word $D'$) — *"the invalid certificate that looks clean"* — caught by a **width-sign check**, the
exact failure mode this target's own earlier SOLVER_AUDIT had predicted.

**An owner vacuous check, confessed to the agent.** The owner's first improvement check looked for
field names that do not exist in the artifact, compared **zero rows**, and reported zero failures —
which reads identically to success. It was only caught because "0 rows improved" was inconsistent
with "0 rows failing." **This is the same disease as the agent's own bug (iii)**, where the
$+\infty$ branch was re-verified on orbit-**aggregated** masses that never vanish, at a granularity
structurally unable to see the thing under test. The agent caught its own; the owner nearly shipped
his. Inconsistency between two derived counts is a cheap and effective guard — adopted.

**`mm3/` off-diagonal, `paper55`: certified minimum 55 over ALL 5,700 survivor data-triples.**
Decided in 929 s. Histogram owner-verified to sum to exactly 5,700; **minimum 55 attained by exactly
one data-triple — the all-monomial one, the known scheme; ZERO at $\le54$; median 66.** So the low
tail is a single point and it is the published result. The certified no-go extends from one monomial
point to $1{,}891{,}123{,}200$ scheme-instances across the $\sigma$-classes.

**A false alarm the agent raised against itself, caught by the anchor.** Its `total` field held the
per-side sum without the $+14$ transposition gap, so the live $\le54$ trigger fired on side-sums
(min 41). **What exposed it: the all-monomial data-triple must certify to exactly 55 and the code
printed 41 for it.** Nothing else was inconsistent — a stream of sub-54 totals is precisely what a
genuine discovery looks like. **Only a value known in advance could separate discovery from
bookkeeping error.** Second instance today, after `omega`'s V1 gate catching a $17.7\times$
aggregation error. Instruction issued for `sun56`: **assert the anchor at startup, do not check it
afterwards.**

**`oct-rank/` S3 filed as a negative with a quantified obstruction** — best certified LB 13, gap to
14 exactly 1, and **the substitution route proven unable to exceed 13** (exact 181-checks-per-triple).
Owner re-verified the Cayley-Dickson convention ($L_{\bar u}L_u=N(u)I$ exact), commutator rank 8 on
five independent basis triples, and that the flattening caps $3/8/8$ are **structural** for an
$8\times8\times3$ tensor. The agent labelled its own most attractive route down twice —
`COMPUTATIONAL-EVIDENCE` for a 173/173 sweep rather than a universal quantifier, and
`CITED-DEPENDENCY, UNVERIFIED` for the Strassen twin form it had not read — then **falsified the
attractive reading itself** with a $\tau$ control at $n{=}4$. Those labels are what kept a
published-bound improvement claim out of the record.

### `mm3/` audit RE-OPENED and the owner's own closure retracted, 2026-08-30 — miss #13, and the sharpest one

**Sequence, because the shape is the lesson.** The owner narrowed a claim, demanded a measurement
rather than an argument, received one, closed the audit and restored the claim — and the closure was
**wrong**, because the measurement's *reference* had not been validated.

- **A1** compared the frozen map against the agent's **first** derivation of the honest action:
  frozen $6912/0/48$ vs honest $6912/0/48$, survivor sets bit-for-bit identical. Owner accepted it
  and restored *"certified $\ge55$ over $S$."*
- **A2**, the agent's own falsifiable-both-ways control, then **disproved that first derivation** —
  it failed Brent 729 on **10/10** random non-monomial triples, and even disagreed with the frozen
  code on monomials. **A1 had compared two non-automorphisms.**
- **Owner miss #13, and it is worse than the previous twelve** because it happened while enforcing
  the exact rule it broke: the owner had written to this agent hours earlier that *"two
  implementations agreeing on the same defective semantics is not corroboration."* Then accepted
  precisely that. **Demanding a measurement is not sufficient — the measurement's REFERENCE must be
  validated first.**

**The corrected family, and the owner's error in the premise.** The matmul trilinear form is
$F(U,V,W)=\sum U_{ik}V_{kj}W_{ij}=\langle UV,W\rangle=\mathrm{tr}((UV)^{\mathsf T}W)$, **not**
$\mathrm{tr}(UVW)$. Corrected: $u'=P^{-\mathsf T}UQ^{-1}$, $v'=QVR^{-1}$, $w'=PWR^{\mathsf T}$,
invariant for **all** invertible $P,Q,R$. **Owner-verified:** symbolic identity with fully symbolic
matrices, plus **200/200** random integer sextuples invariant for the corrected family against
**0/200** for the naive one.

**The reversal.** At the diagonal, corrected gives $u'=G^{-\mathsf T}UG^{-1}$, $v'=GVG^{-1}$,
$w'=GWG^{\mathsf T}$ — the $v$-side differs from frozen's $GVG^{\mathsf T}$ off the orthogonals.
`paper55`: $48$ under both maps, agree. **`sun56`: $48$ frozen vs $576$ corrected, $528$ of them
non-monomial.** So *"the 48 survivors are exactly the monomials"* is **map-relative** — false for
`sun56` under the true automorphism.

**Standing position.** The **288 monomial decisions stand unconditionally** (valid under any
automorphism, Brent-verified, floor-certified). The **exclusion of the 6912 does not**, and now has
a concrete address: **528 `sun56` non-monomial diagonal orientations, ternary and Brent-valid, real
decompositions inside the ternary alphabet, certified LBs UNKNOWN.** Owner has requested one
explicit witness $(G,U',V',W')$ for independent 729/729 verification before the 528 figure enters
the index as more than an agent report. **No published claim is affected** — the gap is this
repository's coverage, not arXiv:2607.28676.

**Second-order finding recorded:** the record-attack manifest's note that *"each $G$ is a tensor
automorphism with integral inverse"* is **wrong off the orthogonals** — 8/8 random non-monomial
diagonals fail all 729 Brent identities under the frozen map, untested there because the ternarity
gate excluded them first.

**Paper-side scope note, correctly NOT escalated by the agent:** the paper's literal
$(XUY^{-1}, YVZ^{-1}, ZWX^{-1})$ with three independent letters is an automorphism only on the
orthogonal subgroup — the same shape as the code defect — while read **contragrediently** it matches
the true family. Published claims hold under their own conventions, every ladder scheme is ternary,
and published checks only ever ran on monomial/orbit sandwiches. A scope note on reading
"orientation" off-diagonal, **not** an error in anyone's paper, and to be kept at that strength.

### `mm3/` audit opened and closed the same day — SUPERSEDED, retained per rule 5

Full arc, because the shape of it is the point:

1. **Challenge.** `Mm3OffDiag`, reading the frozen code for its own purposes, found
   `gatec_sweep.sandwich` hard-coding `inv_sp` — a monomial-only routine — while
   `census_fmpz_path.py` fed it all 6960.
2. **Owner narrowed the index immediately**, not after deliberation: certified $\ge55$ held only
   over the 288 decided orientations, exclusion of the 6912 marked under audit.
3. **Root cause sharpened by the agent, and it corrected the owner.** The owner's diagnosis
   (*"transpose used where inverse was required"*) was the **symptom**. Matching the implemented
   $u'=X_i^{\mathsf T}UY^{\mathsf T}$ against the tensor-preserving family, the $U$-side forces
   $Q=G^{-1}$ and the $V$-side forces $Q=G^{\mathsf T}$; both hold only if $G^{\mathsf T}=G^{-1}$.
   **Cause: a formula that is a tensor automorphism exactly on the orthogonal subgroup, applied
   outside it.** Owner re-derived this independently from `gatec_sweep.py:97-110` and adopted the
   agent's phrasing verbatim.
4. **Owner refused the rescue as an argument.** The agent showed the intended map coincides with the
   honest action at $(P,Q,R)=(G^{-1},G^{-1},G^{-1})$ and that inversion preserves monomiality, so
   the survivor set should be the 48 either way. Correct reasoning — and per rule 14 still an
   argument, so the narrowing stayed until measured.
5. **Measurement.** Diagonal census re-run over all 6960 under **both** maps, same enumeration
   order: **frozen $6912/0/48$, honest $6912/0/48$, survivor sets identical bit-for-bit (48
   monomial), Brent-fail cell $0$ in all four runs, both decompositions.**
6. **Claim restored, stronger than before.** Certified minimum **55 over the full $S$** reinstated,
   now verified under two distinct maps one of which is provably the correct action — evidence the
   original campaign never had. Audit trail retained inline per rule 5, retitled *"what was audited,
   and why it was not a false alarm."*

**It was not a false alarm.** The code really was outside its domain for 6912 of 6960; the count
survived only because of a structural coincidence specific to the **diagonal** slice. That
coincidence is recorded rather than erased by the happy ending, and it gives **no** comfort
off-diagonal — which is exactly why no frozen code is reusable in the running census.

**Owner-verified design figures** (all re-derived independently): right cosets of the 6 permutation
matrices over the 6960 number **exactly 1160 with every orbit size exactly 6** (free action);
permutations factor out so the predicate depends on **data-pairs only**; $1160^3\times648 =
1{,}011{,}460{,}608{,}000$ **identical** to $6960^3\times3$; and the built-in control
$2^3=8$ data nodes $\to 8^3=512$ triples $\to 512\times648=\mathbf{331{,}776}=48^3\times3$ **exactly**,
forcing the machinery to reproduce an independently-derived figure before reporting anything new.

### Systematic follow-up: permissive-direction audit across all seven targets, 2026-08-30 — no new defect

The `kg` defect (a cap that **understates** a quantity bounded from above, so it can manufacture a
PASS) and the `mm3` defect (a precondition violated at a call site) are two instances of one class:
**errors invisible to arithmetic re-checking.** Both were found by agents reading code for their own
purposes. So the owner audited the class directly rather than waiting for the third instance.

**Scope, stated honestly because the instrument is narrow:** a pattern search over `src/` of all
seven targets for hard caps, clamps, clips and `min(1,\cdot)` forms that could sit inside a
certified path. **Result: no new defect.** Every hit resolved into one of three sound categories:

1. **Mathematically justified in the conservative direction.** `kg/gate_B.py:113-115` bounds
   $\lVert\text{odd}\rVert^2 = 1-\lVert\text{even}\rVert^2 \le 1-\max(0,\lVert\bar a\rVert-\rho)^2$ —
   the upper bound on the target uses the **lower** bound of the intermediate, which is the correct
   direction; the `min(1,\cdot)` on the other end is valid for a unit cubic.
   `kg/gate_c_b5.py:25-26`'s $\min(11/12,C_\lambda/2)$ is a genuine minimum of two bounds.
2. **Float-only, outside the certified path.** `omega/alman25_float.py:179`,
   `omega/vxxz24_float.py:223`.
3. **Documented, flagged, and correctly labelled.** `omega/interval_core.py:193-199` clamps
   $-t\ln t$ to $[0,hi]$ **with the validity argument in the comment** and **raises a flag** so the
   box is recorded support-invalid — explicitly *"we avoid silently clipping."*

**The exemplary case, worth naming as the standard.** `delcap/prc_cert.py:23-28` faces a
$10^{-300}$ clip in the *source paper's own code* that changes the entropy zero-term by an amount
**"we cannot rigorously bound."** Its response: **reuse the paper's exact clipping** so the
reproduction is faithful, then label the output `COMPUTATIONAL-EVIDENCE, {REPRODUCED} numerics, not
a certificate`, and repeat the caveat in `gate_ab_run.py:12-13` and in the row's evidence label.
That is the right handling of an irreducibly non-rigorous input: reproduce it honestly, label it
down, and say why — rather than either inventing a bound or quietly promoting the number.

**What this audit does NOT establish**, since a clean negative is worth exactly its instrument: the
pattern search covers one syntactic class in `src/` only. It would not have caught either defect
actually found today — `mm3`'s was a **call-site precondition**, `kg`'s was a **constant inside a
bound derived from a paraphrase**. Neither is a clamp. The durable practice remains the one from the
`kg` entry: **re-derive an inherited bound's worst case from the primary source**, rather than
reproducing its outputs or grepping for its shape.

### `kg/` gate B — an unsound bound found, and two of this repo's own PASSes retracted, 2026-08-30

Agent `KgGateBClose` did not close the three OPEN bands. It found something better: the
predecessor's envelope was **unsound in the permissive direction**.

**Defect 1.** `j_integrand_bound` capped $|e|\le1$. At the paper's own D.4 reproducing kernel
$p_0=(3-s^2)/\sqrt6$ (`paper_full.txt` line 2036), $|e(0)|=\sqrt6/2=\sqrt{3/2}\approx1.2247$.
**Owner-verified exactly:** $a_0=\sqrt6/3$, $a_2=-\sqrt3/3$, $a_0-a_2/\sqrt2=\sqrt6/2$ identically.
**Defect 2.** 3-D circumradius inflation in place of the paper's 2-D even-box envelope (1972-1981).

**Why the direction is the whole story.** Capping $|e|$ too **low** understates $J$, and
understating $J$ makes $J\le d$ appear to hold. **An unsound bound in the permissive direction can
manufacture a PASS** — categorically worse than a loose bound, which can only fail to certify.
Every band certified under that cap became suspect the moment the defect was found.

**Consequences applied to the index:** $[1.30,1.45]$ **stays PASS** — the mandatory anchor was
re-run under the corrected envelope and re-passed at certified margin $+0.1930$, which is precisely
what the anchor requirement is for. $[3.50,4.083]$ **retracted PASS $\to$ PARTIAL** (corner cell
stuck at certified $-2.0$ to $-4.7\times10^{-5}$). $[4.083,6.0]$ **retracted PASS $\to$ UNVERIFIED**
— not re-run, so its PASS rests on the unsound cap; the agent declined to inherit it, writing
*"predecessor PASS remains their own cover's."* New: a $[6.0,12.0]$ razor cell certified at
$+3.675\times10^{-5}$, the tightest mid-band certification in this repo.

**Not a crisis, and established rather than assumed.** Float probes at the stuck corner give
$\max_p J(3.75,p)=0.14467$ vs $d(3.75)=0.15684$ — true spare $+0.0122$, about **600$\times$** the
paper's typical margin. The stuck cells are our envelope excess; **no cell is near refuting Lemma
D.2/D.4**, and zero escalations was correct.

**Three process notes worth keeping.** (i) The agent read the D.3 route **first-hand** with line
numbers rather than from the predecessor's paraphrase — and the paraphrase is exactly where the
2-D even-box structure had been lost, so first-hand reading is what found both defects. (ii) It did
**not** adopt the named panel-subdivision fix mid-campaign, pre-registering it instead — rule 16
held under the pressure that usually breaks it. (iii) It disclosed that Addendum 2 postdates the
`b5 v2` run it governs, **in the addendum itself**. A timing gap declared against oneself is worth
more than a clean-looking timeline.

**Two handoff items from the agent, the second with an owner qualifier attached.**
(i) Its final per-panel weights use **exact $\Phi(p_b)-\Phi(p_a)$ increments** rather than
$\phi(p_b)$ sup weights — both **tighter** (a sup over a panel over-estimates what the exact
increment pays) and **cheaper**, cutting razor-cell panel counts from $>10^4$ to $\approx2\times10^3$.
(ii) It names the razor cell $(0.8,-0.6)\pm10^{-4}$ at $c$-pair $(6.0,6.05)$ as the worst-case
witness for $[4.083,6.0]$, so its certified $+3.675\times10^{-5}$ *"pre-prices the band."*
**Owner qualifier:** "worst case in the band" is a **monotonicity claim in $c$**, and a certified
margin at the band's right endpoint bounds the band only if that monotonicity is itself certified.
Recorded as the **starting point and likely binding cell, not a bound.** This matters precisely
because $[4.083,6.0]$ is the band just retracted to UNVERIFIED: reading the pointer as a bound
would re-certify a band from one cell plus an unproven inequality — reproducing, in a new form, the
very defect this campaign found.

**Repo-level lesson, and it is the second instance today.** Both `mm3` and `kg` had a defect that
**no amount of re-checking the reported numbers could find** — one a precondition violation at a
call site, one an unsound constant inside a bound. Both were found by an agent **reading code or
source it needed for its own purpose**. Owner arithmetic re-derivation is necessary and has caught
real errors today, but it has **zero power against wrong semantics**. The generalisable practice:
every campaign that inherits a predecessor's implementation must re-derive the inherited bound's
**worst case** from the primary source, not reproduce its outputs.

### `omega/` gate C stage 3 — FAILURE TO CERTIFY, reported as such, 2026-08-30 (agent `OmegaSlope`)

**Outcome: the machinery is incomplete, the 21-dimensional nullspace question is UNCHANGED and
OPEN, and no certified number about $\omega$ was produced.** Certified end-to-end coverage: **zero
decisions.** No domain move, no re-centering, no budget request, no patch. Pre-statement committed
(`83e3f09`) before first compute.

**Built and working** — the thing the predecessor named as missing: a certified forward-mode
interval-AD core plus a certified replay of the transcribed VXXZ24 tree carrying **45-slope
bundles through every consumed quantity** (glob, 135 L3 parts, 1377 L2 parts, Lemma-1 $\epsilon$
terms, `part_frac` chains), running end to end at $\approx$110 s per pass for both point and box.

**Failed at value re-aggregation, caught by the anchor gate.** The replay never builds the level-2
`PartLv2 mat_size_contribution`, so `M_low` collapsed $2.09425\to0.14985$ and the point pass
returned $\Omega_{\rm raw}=41.9962$ against the frozen raw $2.3715538358350803$. Every number from
those runs was declared **non-evidence** and entered no claim.

**Owner verification.** (i) The anchor is correctly specified and I checked it against a suspicion
of my own: the agent compared to the **raw** endpoint, not the certified one. `omega/README.md`
lines 171-174 give certified $=$ raw $+\,1.9380875689560958\times10^{-11}$, and raw $+$ defect
reproduces the certified endpoint to $\approx1.2\times10^{-16}$ — the float display granularity the
README flags at line 287. A point re-aggregation must hit the raw value, so the V1 gate is right.
(ii) One phrasing conflation corrected: $14\times$ is the `M_low` ratio ($13.9756$); the endpoint
inflation is $41.9962/2.3715538=17.7083\times$. Two different factors; both figures appear in the
report, so nothing is wrong, but the sentence merges them.

**Why the anchor gate is the load-bearing part.** A $17.7\times$ miss is unmissable. Had the missing
term instead produced a $1.001\times$ error, the slope bundles would have looked healthy, the LP
would have run, and the campaign would have produced a **certified-looking local-optimality
statement resting on a silently wrong aggregation** — the exact failure that has bitten three
targets here. The HiGHS trap was also respected via the u-form, which mattered: the quantity in
play is $\approx1.58\times10^{-7}$ against a default primal tolerance of $10^{-7}$.

**Net movement:** the blocker went from *"new machinery, not built"* to *"one aggregation term,
formula specified, V1 tolerance $5\times10^{-13}$."* Gate C stays PARTIAL; the scoped negative over
the six pre-registered probes stands unchanged; the index is untouched, which is the honest
consequence of producing no certified number.

### `mm3/` frozen-artifact defect found by an agent and confirmed by the owner, 2026-08-30 — a claim recorded hours earlier is now narrowed

Agent `Mm3OffDiag`, reading the frozen record-attack code before designing its own census, found
that `gatec_sweep.inv_sp()` implements the sandwich inverse as the **TRANSPOSE**.

**Owner verification, two measurements:**
1. $G^{\mathsf T}=G^{-1}$ holds for **exactly 48 of 6960**, and that set is **bit-for-bit identical
   to the monomial set**. The reason is structural: an integer matrix satisfies $GG^{\mathsf T}=I$
   iff it is orthogonal, and the integer orthogonal $3\times3$ matrices are precisely the signed
   permutations. So transpose-as-inverse is **exact on the monomials and wrong on all 6912
   non-monomials**.
2. The exclusion cannot be rescued by an inverse-ternarity argument: **4656 of 6960 have both $G$
   and $G^{-1}$ ternary, of which 4608 are NON-monomial.** Worked example:
   $G=[[-1,-1,-1],[-1,-1,0],[-1,0,-1]]$, $G^{-1}=[[1,-1,-1],[-1,0,1],[-1,1,0]]$ — both ternary,
   $G$ non-monomial.

**Consequence, applied to the index immediately.** The certified $\ge55$ **stands over the 288
decided orientations** — those are monomial, so transpose *was* the inverse there and they were
decided under the correct action. The **exclusion of the 6912 is marked UNDER AUDIT**, and the
headline is narrowed from *"certified minimum 55 over $S$"* to *"over the 288 decided
orientations."* If any of the 4608 preserves the alphabet on the real decomposition, $S$ was
smaller than claimed and $\le54$ could hide in the wrongly-excluded region. **No published claim
moves either way** — the record 55 was never beaten here.

**What this says about the process, and it is not flattering to the owner.** The record-attack
result was accepted this session after the owner independently brute-forced **every count in it**
— 11808/6960/48/6912, $41760$, $288$, $145$, $6960^3\cdot3$ — and all were exact. The counts were
right and the **map** was wrong, so no amount of count-checking could have caught it. What caught
it was a sibling agent **reading the frozen code** rather than the frozen numbers, because it
needed the same routine for a different purpose. Verifying arithmetic is not verifying semantics,
and the owner's re-derivation had no power against this defect — rule 14 at the level of an audit
rather than an experiment.

**Also recorded:** the owner's 4608 figure is a **necessary-condition screen**, not a survivor
count — ternarity of $G^{-1}$ as a matrix does not decide whether $G$ preserves the alphabet on the
actual decomposition. The agent was told explicitly not to inherit it as a count.

### Two more gates opened on the idle targets, 2026-08-30 — and the owner falsified his own hypothesis before the agent could act on it

Seven campaigns now run concurrently. Two were opened on the only targets without an agent:

- **`oct-rank/` S3** (`OctRankS3`) — the one open item this target names: is
  $\mathrm{rank}(L_u,L_v,L_w)\ge14$ for independent $u,v,w\in\mathbb R^8$? The chain gives
  $5+13=18$; a $\ge14$ bound gives **19**, improving the published lower bound of arXiv:2608.16649.
  Escalate-before-writing. Anchors required first: $n{=}2\Rightarrow3$, $n{=}4\Rightarrow8$.
  Non-certificate-grade material stays in `scratch/` as `HUMAN-AUDIT-PENDING` per pre-statement C3.
- **`mm3/` off-diagonal** (`Mm3OffDiag`) — the $6960^3\times3=1{,}011{,}460{,}608{,}000$ sandwich
  region named by the record attack as the only remaining place $\le54$ could hide inside the
  ternary alphabet.

**OWNER HYPOTHESIS FALSIFIED BY THE OWNER, BEFORE COMPUTE — miss #12.** The `mm3` brief carried a
hypothesis for making $1.011\times10^{12}$ tractable: that alphabet preservation might force each
of $X,Y,Z$ **independently** monomial, collapsing the space to $48^3\times3=331{,}776$. It was
labelled as a hypothesis to test, and then tested immediately rather than left for the agent:

**300,000 random pairs of non-monomial ternary unimodular matrices, exact integer products:
70,639 — $\mathbf{23.55\%}$ — stay inside $\{-1,0,1\}$, and all 70,639 also have $|\det|=1$.**
Explicit witness where two **non-monomials** multiply to a **monomial**:
$X=[[0,1,0],[1,1,0],[-1,-1,-1]]$, $Y=[[-1,0,1],[0,0,-1],[1,-1,0]]$, $XY=[[0,0,-1],[-1,0,0],[0,1,0]]$,
a signed permutation matrix.

So cancellation between independent factors is **generic**, not exceptional, in exactly the
semigroup the sandwich draws from. The diagonal census's $6960\to48$ collapse is a fact about
$X=Y=Z$ and does not extend. The agent was told within minutes, before designing its census.

**The falsification's own scope, stated so it is not over-trusted either:** the test used plain
products $XY$ as a **proxy**, not the real sandwich action on the decomposition's factor matrices.
It demolishes the *mechanism* the owner invoked; it does **not** determine the true off-diagonal
survivor count, which only the census on the actual action can settle. The agent's brief was
amended to plan for survivors far exceeding $48^3$, with the deliverable restated as *exact
survivor count + certified verdict over a named subset + precisely named unswept remainder*.

**Pattern now at three of the last four misses being attribution or mechanism rather than
arithmetic.** The owner's numbers survive checking; the owner's *reasons* do not. The operational
response adopted today: test every "because X" clause in an owner brief immediately, and prefer
falsifying it before the agent commits runway.

### `mm3/` record attack closed as a certified negative, 2026-08-30 — and it corrected the owner's framing

Gate C's record attack, listed NOT STARTED since the budget wrapped, is now **complete**:
**certified $\ge55$ over a pre-registered named set of 41,760 exactly-decided triples, no
$\le54$ anywhere.** Pre-statement committed before compute. Published record 55 not beaten,
nothing escalated.

**Owner re-verified every count by independent brute force** over all $3^9$ ternary matrices:
11808 invertible / 6960 unimodular / 48 monomial / 6912 non-monomial, $41760=2\cdot6960\cdot3$,
$288=2\cdot48\cdot3$, $145=6960/48$, $6960^3\cdot3=1{,}011{,}460{,}608{,}000$. All exact.

**The agent corrected the owner on attribution, and was right.** The commissioning brief said the
monomial-transfer theorem makes the record attack tractable by collapsing the orientation space.
It does not. The $145\times$ collapse is done by the **ternary alphabet** — 6912 of 6960
non-monomial $G$ leave $\{-1,0,1\}$ under the diagonal action — and is an instance-exact census,
**not a theorem**. The theorem covers exactly the 48 monomials. Owner error; the agent declined to
inherit it. **This is the eleventh owner miss today, and the second of the attribution kind rather
than the hypothesis kind** — the pattern is now specific enough to name: when the owner supplies a
*reason why* something is tractable, that reason is the part most likely to be wrong.

**Agent self-correction under rule 5, owner-confirmed:** *"ternary unimodular group"* is a
misnomer — the 6960 set is **not closed** under multiplication; the owner exhibited a product with
an entry of **3**. The 48 monomials are closed, forming $B_3$, order 48. Only the noun was wrong;
the per-element action is sound.

**Certificates:** 15 UNSAT instances, kissat 4.0.4 DRAT→LRAT accepted by **both** `drat-trim` and
`lrat-check`, plus CaDiCaL 3.0.1 cross-verification; anchors 55/58/56/59/60 reproduced **first** at
729/729 each; two independent arithmetic backends agree on every verdict. All checks over
$\mathbb Z$, so the arXiv:2607.29291 ring hazard is closed by construction rather than by argument.

**Left open and named, not omitted:** off-diagonal $6960^3\times3\approx1.011\times10^{12}$ — the
only remaining hiding place for $\le54$ inside the ternary alphabet.

### Live supervision catch that came back clean, 2026-08-30 — and the dead reduction is now separated by measurement

`DelcapNextGate`'s status line read *"append selftest and main to orbit module."* This target's
own history contains a **falsified** orbit hypothesis — the type-channel symmetry reduction,
killed by decisive numerical disproof because $S_n$-equivariance of the deletion-channel law
fails. An orbit module could unknowingly re-implement it, so the agent was asked one question
before building further: **which group acts on which object, and does the action commute with the
channel?**

**Answer: case (a), and it is correct.** The group is $G=S_q\times C_2$ acting on **symbol
values** — alphabet permutation $\sigma$ plus word reversal $\rho$, applied simultaneously to
input and output words. Deletion acts on **positions**; value relabeling and position deletion
commute, so $A(\sigma x,\sigma y)=A(x,y)$. The agent had already named the falsified reduction
explicitly in the module docstring and pre-statement so it cannot be confused with this one.

**Owner independent verification, exact rationals, deletion channel rebuilt from scratch**
(instrument sanity first: every row sums to exactly 1):

| $(n,q,d)$ | rows sum 1 | $S_q$-equivariant | reversal-equivariant | **control: position swap** |
|---|---|---|---|---|
| $(3,2,1/3)$ | yes | **exact** | **exact** | **BREAKS** |
| $(4,2,1/5)$ | yes | **exact** | **exact** | **BREAKS** |
| $(3,3,2/7)$ | yes | **exact** | **exact** | **BREAKS** |

Exact over all $x$, all $y$, all $\sigma$. **The control is the point:** swapping two input
*positions* breaks equivariance on all three instances, in the same code path that confirms the
value action holds. So the live reduction and the dead one are now separated **by measurement from
both sides** rather than by assertion — value relabeling commutes with deletion, position
permutation does not, because deletion preserves order.

**What earned the check its keep:** the agent had built a **KL-invariance test with a
non-invariant control that the test can detect**, unprompted. A symmetry test with no failing
control passes on a wrong group. That is rule 14 discharged rather than recited, and it is the
same instrument-power question that four separate failures turned on today — twice inside the
owner's own reasoning.

### A consistency check that correctly found nothing to fix, 2026-08-30

Swept every `0-for-N` owner-calibration figure across `README.md`, `RESULTS.md`, `PROGRESS.md`
and the scan doc: eight distinct values (`0-for-3` through `0-for-10`) across twelve locations.
The obvious reading is drift — twelve copies of one number disagreeing.

**It is not drift, and normalizing it would have been the error.** These are a **running tally**
in a **dated, newest-first ledger**, and each figure is correct *as of the entry that carries
it*. Rewriting the historical entries to the current value would destroy the very thing the
ledger exists to record: that the count rose, entry by entry, as each hypothesis died. Exactly
one location carries current-state phrasing (*"now stands at"*), and only that one was updated,
`0-for-6` → `0-for-10` — it had genuinely fallen four behind the ledger.

Recorded so a future consistency pass — mine included — does not "fix" a correct record. **A
count that appears with different values across a chronological ledger is evidence the ledger is
working, not evidence it is broken.** The distinction is whether the number describes *state*
(one owner, must agree) or *history* (many entries, must not be touched). This is the inverse of
the `physics/na-compiler` and `delcap` index drifts closed earlier today, where two copies of a
*state* disagreed and one had to go.

### Index self-audit against the newly-added rules, 2026-08-30 — clean, and one instrument discarded

Held `cs/RESULTS.md` to rules 7 and 15 (name the certified domain) on the phrasing that actually
went wrong today: *locally optimal*, *no feasible*, *exhaustive*, *optimal*, *impossible*,
*cannot be beaten*, *no X exists*.

**Result: 16 high-risk claim sentences, zero violations.** Five needed hand-reading and all five
resolved benign — *"float-optimal point"* and *"an honest enclosure of a float-optimal point
costs…"* are **descriptive** (optimal in the float solver's sense, which is the honest
qualification, not a universality claim); *"floor-impossibility results"* is the **name** of a
result type whose domain lives in the section referenced on the same line; *"exhaustive
line-support minimization"* names its own domain. The remaining eleven carry an explicit domain
marker in the same sentence.

**The most reassuring hit is the owner's own retraction**, surfaced by the audit sitting inline
where rule 5 puts it: *"The agent's original wording was the weaker 'certified-locally-optimal …
within every feasible single-block family **tried**'; I sharpened it into the universal form,
which is the stronger and unsupported claim."* The audit found the error already recorded against
the person who made it.

**One instrument discarded rather than reported.** A first, broader pass flagged **45**
universal-language sentences with no domain marker. Hand-reading the first dozen showed it was
near-worthless: *"all seven targets"*, *"to repo knowledge"*, *"bit-identical at all five"*
(preceded by *"at MID precisions 300/200/128/96/64"*), *"not constant on that kernel"* — every one
carries its domain, and the regex simply missed *"on that"* versus *"on the"*. **That count is not
reported as a finding**, because reporting it would have been precisely the fabricated-drift
failure recorded one entry above. Two instruments, one signal: the index is clean on the phrasing
that has actually caused errors here.

### Forward coverage sweep, 2026-08-30 — ceiling moved, no target affected

The morning audit swept backwards and found the ECCC report-number ceiling defect. This sweep
went forwards. **arXiv `cs.CC`/`cs.IT`: nothing new and the absence is verified** — newest
announced submission 2026-08-27, and 08-29/08-30 are a weekend with no announcements; the max
`published` date was printed before concluding, per rule 17b. Next window Monday 2026-08-31.
**ECCC: one new report past our ceiling** — `TR26-162` (Xiao, *Weighted Bipartite Matching is in*
$\mathrm{Mod}_p\mathsf{L}$), paper-dated 29 Aug 14:59, **publication-dated 30 Aug 18:05**, so it
landed after both the scan and the morning audit. Triaged against all seven targets: **no
collision, no gate affected, no target opened.** Ceiling `TR26-161` → `TR26-162`.

**And the coverage coordinate is still not precise enough.** ECCC carries a paper date *and* a
publication date, differing by over a day here. This morning's fix ("track by date, not report
number") therefore leaves a second silent ceiling unless the sweep names **which** date it
bounded by: `TR26-162` is invisible to a paper-date sweep run before 18:05 today and visible to
a publication-date one. Two ceilings found in one day on the same source, both by looking rather
than assuming.

### Cross-repo SSOT audit closed, 2026-08-30 — and the owner's own instrument failed twice

All three sibling repositories verified SSOT-clean:

| repo | targets | index coverage | verdict |
|---|---|---|---|
| `cs/` | 7 | 7 rows | clean |
| `physics/` | 6 | 6 rows | clean (the `na-compiler` gap closed earlier today) |
| `math/` | 13 | `uc/` + `LIU_H1/` only, **declared** | clean — scope stated, not implied |

**`math/` is not a drift and the owner's framing of it was wrong.** `math/README.md:3-12`
explicitly declares `RESULTS.md` as the authoritative index *for the union-closed programme
only*, names the sibling difference, and instructs readers not to cite it outside
`uc/`/`LIU_H1/`. Path-anchored measurement confirms the declaration is **exactly** honored:
214 `uc/` references, 14 `LIU_H1/` references, and **zero** references to any of the other
eleven targets. Stating a narrower scope is a better resolution than manufacturing one row
per target, which is what the owner had proposed — a document that grew as one project's
own index should say so rather than pretend to a coverage it does not have.

**Two instrument failures in a single owner check, both caught before anything was recorded:**
1. **Wrong standard.** The check tested `math/` against one-row-per-target — a convention
   `math/` explicitly rejects in writing. It reported 13 targets "MISSING from index" for a
   document that disclaims them by design. An audit that flags a *declared* difference as a
   defect is measuring the wrong thing.
2. **Unanchored pattern.** Target-mention counts used bare substrings, so `ns` scored **140**
   and `dist` **11** — matching `ns` inside *constraints*, *functions*, *means*, and `dist`
   inside *distance*, *distinct*, *distribution*. Path-anchored, both are **0**. The owner was
   one step from recording "`ns` has more mentions than `uc`, so the declared scope is
   violated" — a fabricated drift produced entirely by the measurement.

**The sharper point, and the reason this is recorded rather than quietly fixed:** the
session-start observation that opened this thread — *"62 `uc/` references, zero for the other
nine targets"* — was made with **the same unanchored instrument**. It was right in direction
by luck and wrong in detail: `LIU_H1/` is not zero (14 references, and it is an in-scope
target), and the repository has 13 targets, not nine. A conclusion that happens to be correct
does not retroactively validate the method that produced it.

- **Discipline rule added (README 17b):** a count is not evidence until its pattern is anchored.

## 2026-08-30 — owner post-scan audit: sources, a citation hazard, a scan gap; six gates launched

Owner-level work while the gate agents run. Nothing here rests on an agent report.

- **Source-version sweep (arXiv API): no source superseded.** `2608.16884`, `2607.28676`,
  `2608.16649`, `2305.07156` all still **v1**. `2608.11158` ($K_G$) is at **v2, updated
  2026-08-12**, which precedes the 2026-08-29 owner read — so we read v2 and the headline is
  not stale.
- **Citation hazard found in the headline's source, and locked conservatively.** The v2
  *abstract* of `2608.11158` gives the upper bound only as the closed form
  $\pi/(2\log(1+\sqrt2))-10^{-4}$. Arb at 300 bits: that is $1.7821139781913691118$, i.e.
  **weaker than the paper's own certificate** $1.7818666069360661$. Full ordering verified:
  ours $1.78184132<$ certificate $1.78186661<$ abstract $1.78211398<$ Krivine $1.78221398$.
  So our margin is $2.5283056\times10^{-5}$ against the certificate but
  $2.7265431\times10^{-4}$ against the abstract — a **$10.8\times$ free inflation** by citing
  a different line of the same paper. `RESULTS.md` cites the certificate and now says
  explicitly that this is deliberate. Comparing against a source's weakest published
  statement rather than its best is how a real result becomes an overclaim.
- **Scan gap: one paper missed, triaged, no collision.** **arXiv:2608.27434v1**
  (2026-08-27, Kassabov–Landsberg–Souza–Speegle, centroids of matmul tensors) appears
  nowhere in the scan or its raw scout output. It refutes a conjecture that tensors with very
  large centroids do not exist, and builds better laser-method tensors from old ones:
  Strassen $q{=}4$ $2.48289\to2.46016$, $q{=}5$ $2.47849\to2.46710$, Schönhage
  $2.548\to2.522$, a $q{=}6$ case at $2.416$ and $2.39$. All far above the record; the paper
  itself states the record as $\omega<2.371177$. **No threat to `omega/` or `mm3/`.** Cause
  of the miss: the scouts covered cs.CC 2026-08 as a listing and this is a final-days entry —
  the month-boundary resolution was the weak point, not a keyword set. Forwarded to both
  affected agents.
- **`omega/` record-artifact audit re-run — still negative.**
  `google-deepmind/alphaevolve_results` holds only the old notebook, last pushed
  **2026-01-05**; GitHub repo searches for `matrix multiplication exponent alphaevolve`,
  `2.371177`, and `matrix multiplication exponent verification laser method` all return
  `total_count = 0`. Still a statement about artifact availability on a date, not about the
  bound. **Consequence:** the record has no reported optimum to box, so gate C's certified
  local optimality must center on a rung with published parameters and be scoped under rule 7
  as a statement about that rung. Sent to `OmegaGateC` before it pre-registered.
- **ECCC coverage defect found in our own scan: a numeric ceiling is not a date ceiling.** The
  scan states ECCC coverage as `TR26-001..154`. ECCC report numbers are **not chronological** —
  across TR26-142..161, **11 of 19 adjacent number pairs are date-inverted**; TR26-154 is dated
  28 Aug while TR26-161 is dated **12 Aug**, 17 days *before* the scan. Seven reports dated
  12–28 August (TR26-155..161) were excluded purely by number. Triaged, none colliding with an
  open target; three are strong future candidates on this repo's own criterion of finite
  checkability — TR26-159 *Sorting from Counterexamples* (Alon–Moran–Moran; exact query
  complexity computable by exhaustive game-tree search at small $n$), TR26-156 *3-Restricted
  Matching Vector Families* (integer-verifiable objects in $\mathbb Z_m^n$ driving
  constant-query LDCs), TR26-157 *Tight Cycle-Cover Inequality for SCS*
  (Chukhin–Kulikov–Mihajlin–Smal; tightness refutable by one small counterexample).
  **TR26-160** (Golowich–Tamo–Zhu, transversal non-Clifford gates from cup products) belongs to
  `../physics/` — forwarded to its ledger, adjacent to `msd/` and `qldpc-dec/`. **No target
  opened here.**
- **Methodological rule adopted, binding on future scans: state source coverage by DATE against
  the dated listing, never by report or identifier number.** Both misses this session are the
  same error in different clothing — a coverage claim expressed in a coordinate the source does
  not order by (report number for ECCC, month boundary for arXiv).
- **Applied the new rule to the scan's own arXiv claim, and it found a second, worse miss.**
  Re-derived the August cs.CC set from the API by submission date: **136** entries against the
  scan's stated 140 from the monthly listing (announcement-month vs submission-date bookkeeping,
  not a gap). Ten matched a target keyword; eight already in-repo; the two that were not are
  `2608.27434` above and — the serious one — **`arXiv:2608.14817`, *"The König constant is
  one"* (2026-08-14), with zero mentions anywhere in `cs/`.** It proves
  $\mathfrak K_{\mathrm K}=1$ for the König bilinear form, states that an elementary Fourier
  argument gives $\le1$ and **excludes equality for every finite dimension**, and gives a
  negative answer to BMMN's high-dimensional question on determining $K_G$ through *alternating
  Krivine rounding schemes*. `kg/` gate C(a) enlarges a **truncated** Krivine-type scheme and
  admits "a certified cap showing the enlarged family cannot beat it" as a valid outcome — and a
  sup equal to $1$ with equality excluded in every finite dimension is exactly that shape of
  cap. Whether it binds *this* family is a question of fact, put to `KgGateC-2` as
  yes/no/cannot-tell with required line citations, **to be answered before the 200k-box budget
  is spent**. Two failure modes rejected in advance: closing a gate on an abstract, and
  reshaping a box to agree with an unread paper — numerology with a citation attached.
- **The citation hazard is a property of the whole source family, not one sentence.** The
  companion case study `2608.11195` restates the same conservative
  $\pi/(2\log(1+\sqrt2))-10^{-4}$. Both the v2 abstract and the companion publish the weaker
  $10^{-4}$ form while the shipped certificate is sharper, so citing the certificate is the only
  comparison that survives reading the full family.
- **Owner-read the missed paper instead of delegating it, and its bibliography exposed four
  2026 $K_G$ preprints the repo did not know.** `Hei26b` **arXiv:2606.00247** (Heilman, *An
  Upper Bound on Grothendieck's Constant*) — bounds $K_G$ by thresholding a **degree-five
  Hermite polynomial** on a random plane, resolving a BMMN 2011 conjecture, with a corollary
  $\text{Kri}-10^{-217}$ at degree three and **"a rigorous computer-assisted proof that
  $K_G<\text{Kri}-10^{-5}$ using interval arithmetic"**. `Hei26a` **2603.22616** — $K_G\ge
  c+10^{-26}$. `JM26` **2603.30039** (Jones–Malavolta) — $K_G\ge K_{\rm DR}+10^{-12}$, because
  every near-extremizer of the Davie–Reeds problem carries $\Omega(1)$ weight on its
  **degree-3 Hermite coefficients** so a small **cubic perturbation** raises the integrality
  gap. `LSX+26` **2606.03991v3** — $\text{Kri}-10^{-5}$, by the **same group** as our source
  (`2608.11158` is `[SLX+26]` in that bibliography).
- **Headline survives and is the strongest bound now known.** Arb at 400 bits, weakest first:
  Krivine $1.782213978191369112$; $-10^{-500}$ and $-10^{-217}$ both numerically invisible at
  this scale; Heilman-interval and LSX+26 $1.782203978191369112$; `2608.11158` abstract
  $1.782113978191369112$; its certificate $1.781866606936066100$; **ours
  $1.781841323880437288$** — $3.62654311\times10^{-4}$ below the best previously published
  *rigorous computer-assisted* bound. Lower side: $6\pi/11$ dominates both 2026 $K_{\rm DR}$
  improvements by $0.0366393187$.
- **A real methodological collision, and `kg/` gate C is HELD because of it.** Heilman reaches
  his bound by interval arithmetic on Hermite thresholding at **degrees three and five** — the
  same instrument and the same degrees as gate C(a) (enlarge the Hermite family, certify by
  interval B&B) and gate C(b) (**the degree-five analogue**). JM26's cubic-perturbation
  mechanism sits in C(b)'s territory too. Neither the `kg/` contract nor its pre-statement
  cites any of them. The agent must answer in writing — is Heilman's degree-five family the one
  we truncate, does his computation already cover our box, is C(b) duplicated/adjacent/independent
  — **before** the 200k-box budget is spent. A clean "already done by Heilman, here is the line"
  is an acceptable and valuable outcome; shifting the box to an untouched corner to manufacture
  novelty is not.
- **Owner hypothesis retracted, rule 5 — now 0-for-5.** On the abstract alone the owner told the
  agent that $\mathfrak K_{\mathrm K}=1$ "may be exactly the shape of a cap" on truncated
  families. **Wrong.** The paper's §1.3 is *"Application to improving the Grothendieck
  constant"* and Corollary 1.3 gives a randomized Krivine scheme with $K_G<\pi/(2\log(1+\sqrt2))$
  — an improvement. The direction runs the other way: a *larger* König value gives a *better*
  bound, which is how BMMN's planar counterexample became the first strict improvement over
  Krivine. Decisively, [NR14] Naor–Regev, *Krivine schemes are optimal*, showed
  finite-dimensional oblivious Krivine schemes approximate $K_G$ **arbitrarily well** — truncated
  families are not capped at all, and that theorem is the licence for gate C(a). Corrected to the
  agent in the same message that carried the retraction.
- **A four-decade retrieval failure looks like a year error.** `kg/README.md` carries
  $K_{\rm DR}=1.6769\ldots$ `[REPORTED]` with Reeds' `bound2.dvi` dead and "Davie 1984
  unobtainable". Two independent 2026 sources now print **16 digits**
  ($1.676956674215576\ldots$), and `2608.14817` cites it as **[Dav85]** A. M. Davie, *Matrix
  norms related to Grothendieck's inequality*, in *Banach Spaces (Columbia, Mo., **1984**)*, LNM
  **1166**, pp. 22–26, Springer, **1985** — 1984 is the *conference* year, not the publication
  year. Owner-verified in Arb at 300 bits: $\pi/2<K_{\rm DR}<6\pi/11<\text{ours}<\text{Krivine}
  <\sinh(\pi/2)$. Stays `CITED-DEPENDENCY` regardless of retrieval.
- **Generalized the bibliography-mining technique to the other targets; two more unknown papers,
  one load-bearing.** `2608.16649` (oct-rank) cites nothing unknown. `2607.28676` (mm3) cites
  one unknown, `2008.03759` (Karstadt–Schwartz bilinear sparsification — leading coefficients,
  not operation counts, so irrelevant to gate C). A keyword sweep then found **`2607.29291`**
  (2026-07-31, *"SAT Certificates for the Matrix-Multiplication Challenges over F2: All Ten
  `Expected-UNSAT` Instances Are Satisfiable…"*) and **`2605.21738`** (2026-05-20, *"Asymptotic
  Rank Speedup Theorems, Revisited"*), both with zero prior mentions in `cs/`.
- **`2607.29291` both validates `mm3/` gate C and opens a bounded audit of gate B.** It builds
  witnesses using *"the $\mathrm{GL}(3,2)^3$ isotropy action, cyclic trace symmetry, and perfect
  matching of transformed summands to constrained slots"* with *"zero residual in all 729 Brent
  equations"* — i.e. the sandwiching action is a published constructive tool and 729/729 Brent
  is the accepted check, both already required of the agent. The matching step is
  canonicalization: group action permutes the 23 summands, so an order-sensitive counter would
  report relabelling as landscape (flagged as an anchor test). **The warning:** their principal
  finding is that ten formulas the Heule–Kauers–Seidl benchmark expected UNSAT are
  **satisfiable**, because the CNFs *"require selected incidences but do not forbid additional
  type-3 incidences"* — an encoding-completeness defect where every solver run was correct and
  every conclusion wrong. `mm3/` gate B's "three independent ways" (DFS, ILP, SAT+DRAT/LRAT) are
  **three checks of one modelling premise, not three checks of the model** — the same principle
  as the `kg/` finding that audits do not transfer across a defect in the layer they exercise.
  Ordered a bounded written audit (per CNF: which clauses *forbid* vs merely *require*, plus the
  both-directions correspondence between assignments and circuits), escalation not a quiet patch
  if an under-constraint appears. **Not a claim gate B is wrong.** No collision with the $\le54$
  target — their object is rank-23 over $\mathbb F_2$ with a type-3 condition and counts no
  additions; imported $\mathbb F_2$ schemes need their own 729/729 check over $\mathbb Z$.
- **`2605.21738` is a scoping constraint for `omega/`, not a collision.** It generalizes the
  Coppersmith–Winograd/Strassen speedup theorems via "Strassen calculus", putting the asymptotic
  rank of $\mathrm{cw}_2$ below $3.931$ against $\underline{\mathrm R}=4$, and bounding any
  $d\times d\times d$ tensor below $d^{2\omega/3}$. Together with `2608.27434` that is **two
  live 2026 lines that change the construction**, whereas gate C certifies optimality *within* a
  fixed construction — so its rule-7 sentence must name the rung, parameter vector, box radii,
  the unsearched record, and the exclusion of differing constructions.
- **`kg/` gate C duplication finding delivered, and gate C is CLEARED.** `KgGateC-2` owner-read
  both colliding papers in full and answered the three ordered questions with technical
  discriminants. **(a) Not the same family:** Heilman `2606.00247` thresholds a one-plane
  $f(x)=\mathrm{sgn}(x_2-\eta h_d(x_1))$ expanded perturbatively about $q=\eta^2\to0$ around
  $\arcsin$; ours is the two-coordinate $\rho$-correlated scheme with
  $\mathrm{sgn}(\pm\vartheta\psi_3)$ transverse partitions, expanded about $H$'s own series and
  closed through a different sufficiency chain (his $A_q(\gamma)=1$ at exact $\gamma$ versus our
  $\gamma+\Delta_H<b_1$ head). **(b) He does not cover our box:** perturbative bound evaluated at
  $q=e^{-450}$, i.e. $\eta\sim10^{-195}$, plus **one** exact-$\eta$ interval certificate at
  $\eta=0.04249900400783211$, $d=3$ (Sage, 54 h Rouché) — nothing at non-small $\eta$, nothing
  multi-coordinate, no optimizer search. **(c) He does not answer C(b):** his degree five is a
  *threshold polynomial inside $f$*, not a coefficient inequality in the inverse-majorant chain.
  JM26 confirms from the SDP side — $\Pi_1-\lambda I-\varepsilon\Pi_3$ with
  $\mathbb E(\Pi_3f)(\Pi_3g)\ge0.046$, no inverse-majorant chain, no $\gamma$, no $b_5$. **C(b)
  is genuinely open.**
- **A correction that cuts against us, made for consistency.** The agent's full-text read gives
  Heilman's **Theorem 1.9** as $\text{Kri}-1.013\times10^{-5}$ where his *abstract* says
  $10^{-5}$; Arb at 400 bits confirms $1.782203848191369112$ versus $1.782203978191369112$, so his
  theorem is $1.3\times10^{-7}$ **stronger than his own abstract**. `RESULTS.md` now cites the
  theorem. The reason is symmetry of discipline: this session insisted on comparing against
  `2608.11158`'s shipped certificate rather than its weaker $-10^{-4}$ abstract, and quoting a
  competitor's rounded abstract while quoting our source's sharp certificate would silently bias
  every comparison in our favour. Our stated margin below the best published rigorous
  computer-assisted bound drops from $3.62654311\times10^{-4}$ to $3.62524311\times10^{-4}$ —
  $1.3\times10^{-7}$ of claimed margin surrendered for a comparison that survives scrutiny (**corrected 2026-08-30: first written as $1.3\times10^{-9}$, a $100\times$ understatement of what we gave up — the two bounds differ by exactly $1.30\times10^{-7}$, so the margin does too**). Ours
  remains the strongest known bound and the headline is untouched.
- **Literature audit completed across all six targets, with per-target verdicts.** `kg/` — five
  unknown papers, gate C **held then cleared** on the finding above. `mm3/` — two unknown, no
  collision on $\le54$, **gate B encoding audit ordered**. `omega/` — two unknown, no collision,
  both bind gate C's rule-7 scoping because they change *constructions* rather than parameters.
  `delcap/` — no collision; all four gate-C primaries are the current frontier; one adjacent find
  (`2604.07234`) carries a **category trap** flagged to the agent, since its object is the rate
  of *uniformly-random codes*, not $\mathcal C(\mathrm{BDC}_d)$, and its own closing sentence
  concedes that its tighter bracket is about a different quantity. `oct-rank/` — no collision;
  the $18\le\mathrm R_{\mathbb R}(T_{\mathbb O})\le25$ window is uncontested, the other
  octonion-and-rank hits being differential geometry. `rs-pe3d/` — **no collision, verified**:
  arXiv `abs:"product expansion"` returns ten hits and every one is physics (operator product
  expansion), `abs:"private PCP"` returns zero, and no other 2026 ECCC report mentions product
  expansion, so TR26-150 is effectively unique. Told that agent the standard of care goes **up**,
  not down: with no competing literature to catch an error, its pre-registration and rule-7
  boundary are the only safeguards.
- **Both coverage claims now verified by date; one limit remains.** `cs.IT` by submission date:
  **402** for 2026-07 and **351** for 2026-08, **753** total against the claimed **776** from
  listings, a $-3.0\%$ offset that matches `cs.CC`'s $136$-vs-$140$ ($-2.9\%$) almost exactly. A
  consistent proportional offset in the same direction across two independent categories is the
  signature of announcement-month vs submission-date bookkeeping, not of missing entries, so both
  claims stand as restated by date. `delcap/` was additionally checked by targeted phrase search,
  which is the *complete* instrument there — any paper bounding deletion-channel capacity contains
  the phrase — recovering all four gate-C primaries plus `2607.03989` and `2604.07234`. **The one
  unclosed gap:** TR26-150's reference list could not be mined from the ECCC abstract page, so the
  `rs-pe3d/` verdict rests on external searches plus its agent's direct PDF read.
- **The lesson worth keeping.** The original scan missed one paper; reading it rather than filing
  it exposed **five** unknown works, one of which uses this repo's own instrument — interval
  arithmetic on Hermite thresholding — at exactly the degrees `kg/` gate C was about to sweep.
  The miss was one paper; the gap was five.
- **Six un-started gates launched in parallel**, disjoint folder ownership: `omega` gate C
  (certified local optimality), `mm3` gate C (orientation sweep beyond the order-3 orbit —
  the predecessor swept only $\{Id,\sigma,\sigma^2\}$ giving $55/57/57$), `delcap` gate C
  (finite-length rows + $q$-ary sandwich), `kg` gate C (septic extension + degree-5 analogue),
  `oct-rank` gate C (the $n=8$ peeling probe, previously "nothing filed"), `rs-pe3d` gate B
  (the open balance/imbalance sweep). All five that yielded placeholder JSON in under 40 s
  were re-driven with explicit first-call lists; **a placeholder status object is not a
  deliverable** and was rejected as such.

### First gate results in: two owner catches, and one owner overclaim

- **`omega/` gate C PARTIAL on the VXXZ24 rung** — see the retraction subsection below; the wider local-optimality claim written here was withdrawn the same day (campaign
  `2026-08-30T11:24:41Z_17bdaeb6_099b1224d603`). **Part 1 settles the $2.026\times10^{-6}$ gap by
  measurement: it is NOT rounding.** The same interval code path at MID $300/200/128/96/64$ bits
  gives a **bit-identical** certified endpoint at all five precisions, gap fixed at
  $+2.0258544615181506\times10^{-6}$, total rounding contribution bounded by
  $8.3\times10^{-19}$ of $\omega$ even at 64 bits. A gap invariant across four octaves of
  precision is not an arithmetic gap. Attribution: the shipped multipliers leave an ln-space
  Lemma-1 constraint-equality residual of $3.60\times10^{-11}$, **inside VXXZ24's own stated
  $1.1\times10^{-9}$ tolerance** — so the release meets its own specification and the gap is the
  price of a stricter standard. **Stricter standard, not defect**; correctly no escalation.
- **Owner catch: I refused to record three "improvements" until feasibility was answered, and all
  three were infeasible.** Fourteen signed boxes were screened over
  $B=\{p^*+\delta:\lvert\delta_i\rvert\le10^{-7}\}^{6759}$; none certified below the published
  $2.37155181$ (no record claim), but three appeared to certify below the $B_0$ anchor, the best
  closing $40.23\%$ of the gap (owner-recomputed at 300 bits). Because the agent had itself listed
  "feasible-repaired dist-simplex probes" as open, I asked the deciding question before anything
  was written. **Answer: all three INFEASIBLE in exact dyadic arithmetic** — `region_prop+` breaks
  $\sum r_p=1$ by $+3.000000000086\times10^{-7}$, `glob_dist_0+/1+/2+` break their dist-simplex
  sums by $+7.0\times10^{-7}$ and the marginal equality rows. Retracted by the agent inline across
  three campaign files with checksums resealed; certified rung-2 stays $2.3715538358544617$.
  Had the phrase "certify strictly below $B_0$" reached the index unqualified it would have read
  as an improved bound on $\omega$.
- **What replaced it is stronger: an exact no-go.** The $r_p$ block is pinned by **three** rows —
  $\sum r_p=1$, $Y-Z=0$, $X-Y=0$ — forcing a uniform exact triple, so asymmetric offsets break an
  asymmetry row and uniform offsets break the simplex row. Exactly-feasible
  $\lvert\Delta r\rvert\approx5.5\times10^{-17}$, ten orders below the $1.94\times10^{-11}$
  absorption. **No exactly-feasible `rp` move exists at radius $10^{-7}$** — the retracted
  $-8.15\times10^{-7}$ is unrealizable in principle, not merely unrealized. The agent also
  corrected its *own* bad filter mid-audit (a "no $-1$ rows" scan that missed rows 85/86) and
  recorded the correction rather than the conclusion it had supported.
- **A latent defect in a validation accessor — provenance moves, numbers do not.**
  `ParamManager._lin_beq` is **replaced** rather than extended on every `add_lincon_eq` call, so
  `pm.linear_violations()` compares all $3174$ rows to only the **last** $b$. Both frozen rungs
  quoted "linear violations $0.0$" through it. Those claims are **RETAINED** — an independent
  exact-dyadic per-row replay at $p^*$ also returns $0.0$ — but now rest on the replay, not the
  accessor. Fix approved in `src/` only; frozen campaign scripts stay immutable.
- **Third instance this session of one failure shape, now a standing rule.** A validation that
  shares a component with the thing it validates cannot see a defect in that component: `kg/`'s
  finite-difference audit checked a wrong function against the analytic derivative of *that same
  wrong function*; `mm3/`'s three-way floor agreement is three solvers over **one** modelling
  premise; `omega/`'s accessor compared every row to a single $b$. **Independent instrument, or it
  is not an audit.**
- **Owner catch in `mceliece/`: a hole in the ladder at $m=9$.** The agent reported "$m\le11$ is
  ALREADY COMPLETE". On its own enumeration — $(6,64,\{3,4,5\})$, $(7,128,\{3,6\})$,
  $(8,256,\{3,5\})$, $(10,1024,\{24,40\})$, $(11,2048,48)$ — **$m=9$ appears nowhere**, and
  $m=10\,t{=}32$, $m=11\,t\in\{32,64\}$ are on its own not-started list. The universal form is
  withdrawn per rule 5; the claim of record is **ten named instances across
  $m\in\{6,7,8,10,11\}$ at the named $t$**, which is how `RESULTS.md` already words it. Confirmed
  the overclaim never reached disk. **Root cause, agent-reported: the plan followed the README's
  ladder rows, which start at $m=10$, so the contiguous intermediate was never added.**
  Re-prioritized $m=9$ ($n=512$, $t\in\{4,8,12\}$, $k=n-mt$ positive with margin) **ahead of**
  $m=12$: it is one step above the completed $m=8$, cheap, and removes the one visible gap —
  a contiguous $m\in\{6,\dots,11\}$ ladder is a materially stronger object than a holed ladder
  plus one extension, because a skipped $m$ is the first thing a sceptic asks about. $m=12$
  continues under a 2 h compute cap with a rule-7-scoped freeze as fallback. Next frozen artifact
  must carry a two-part boundary table: $(m,n,t)$ actually run with seeds, and a per-$m$ line
  naming every $(n,t)$ **not** run.
- **$m=12$ stopped at the cap and frozen unreached — the cap held.** Compute since $m=9$ froze:
  **1 h 54 m against a 2 h cap**, stopped without requesting an extension. Nothing for $m=12$ is on
  disk beyond the $m\le11$ records, no census row completed, and its honest status is *build reaches
  $1478$ s, no verdict*. Frozen as `COMPUTATIONAL-EVIDENCE-in-progress` with the exact stage reached
  and **explicitly not part of the certified ladder**, which remains **13 instances, $m$ contiguous
  $6$–$11$**.
- **The agent's own bottleneck prediction was wrong, and it said so.** It had named the numpy
  fancy-indexing `array_subscript` line in `Instance` init as the sole remaining cost. It was not
  that: the cost is the $\alpha$-guard evaluating all $k=2912$ coordinates at all $n=3488$ points
  ($\approx34$G gathers). Reported as a miss rather than quietly redirected — which is the only
  reason the new diagnosis is trustworthy.
- **But the optimization it *did* land is a structural result, not just a speedup.** The
  family-restricted block factors exactly as $\Psi=C\!\cdot\!J$ with $J$ carrying only
  $2R+1=9$ jet rows, so the block reduces to $(C\!\cdot\!T)^{\!\top}$ of at most 9 rows instead of a
  $k\times k$ contraction. Regression-tested against the naive construction on 4 frozen instances:
  ranks and stacked $N_{\rm fam}$ identical, speedup $7\times$ at $m=8$ and growing in $k$.
- **A proposed fix that is sound but changes the EVIDENCE LABEL — pre-registration ordered before
  implementation.** Apon's identity is $\mathbb F_2$-**linear** in the coordinate index (the LHS is
  linear, and $G^2(\sum c_jf_j)^2=\sum c_jG^2f_j^2$ for $c_j\in\mathbb F_2$), so the good set is a
  **subspace**; if proper, its index is $\ge2$ and 16 seeded random $\mathbb F_2$-combinations miss a
  failure with probability $\le2^{-16}=1.526\times10^{-5}$ (owner-verified). The saving is real:
  $k\cdot n=10{,}157{,}056$ coordinate evaluations deterministically versus $16\cdot3488=55{,}808$,
  i.e. $182\times$ fewer sweeps. **But each probe being exact does not make the conclusion exact.**
  $m\le11$ used the full deterministic guard and is exact; $m=12$ under a 16-probe guard would be
  Monte Carlo. A ladder carrying **one label across two standards of evidence** is precisely the
  failure this repo exists to prevent, and it would be invisible to a reader. Ordered before any
  code: commit probe count and seeding in advance; label $m=12$ rows `COMPUTATIONAL-EVIDENCE` with
  the one-sided error as a first-class number and **never** list them under one heading with the
  $m\le11$ rows; run the deterministic guard on at least one reduced $m=12$ instance to anchor the
  probabilistic path against the exact one (worth more than additional probes); and leave the
  $m\le11$ label untouched. Declining the second evidence standard entirely and freezing $m=12$ as
  unreached was stated as an equally acceptable answer.
- **ePrint sweep: NOT a clean negative — and the revision landed exactly in the coverage gap the
  owner had flagged that morning.** The owner's source-version audit covered the five arXiv sources
  via the API but explicitly excluded `mceliece/`'s IACR ePrint primaries, noting ePrint "has no
  equivalent dated API I could query the same way". A real in-place revision then appeared there.
  **Vedenev 2026/1747 was revised on 2026-08-30, after our copy** — detected by **byte comparison,
  not by a listing date**: cached 245,667 bytes md5 `cd044d4c3b7e560eeea8ed43a6f75a10` against
  today's 264,802 bytes md5 `f27d10be655ddca288d127c3dacf892e`, text growing 1633 to 1835 lines. The
  addition is a new **Section 9, "Extended multiplicity shortening and Frobenius alignment"**:
  Definition 6 extends multiplicity shortening to scalars in $\mathbb F_{q^m}$ (and unlike the
  base-field version it **does** distinguish Frobenius branches), Definition 7 + Proposition 9 give
  the dual/parity-check form, Corollary 6 the punctured-GRS containment. The payload: reading the
  Frobenius alignment off $\dim_{\mathbb F_{q^m}}D^{*2}$ of an intersection of two extended
  shortenings, which he states **removes the $m^{c-1}$ alignment factor from Step 3's cost**,
  leaving only the $q^{m(c-2)}$ support enumeration; Remark 11 credits Saarinen 2026/1786 for
  observing the factor may be removable.
- **Correctly scoped: the revision does not touch what the census measures.** Section 8.1 is intact —
  $c^*=2\delta+3$ for binary Goppa still asserted (line 1404), Apon's Theorem 1 still quoted as the
  lower bound, and the v2 rejection criterion (24) plus Conjecture 5 unchanged. So the Step-3
  convention gate A adopted remains the disputed object and **the 13-instance ladder stands as
  measured**. What v3 changes is cost framing elsewhere in the attack, which this gate never priced.
  Nothing folded in.
- **A fifth paper exists in the wave, but it is not a fifth hold-out paper.** ePrint **2026/1778**,
  Saarinen, *"HOVER: Higher-Order Vanishing Endomorphism Recovery"* (updated 2026-08-24), attacks via
  Hemmert–Wiemers higher-order vanishing — a **different distinguisher family** from the GIJS
  hold-out route — so the hold-out dispute proper is still four-paper (1630, 1747, 1786, 1810). It
  matters regardless: it reports **end-to-end key recovery on five TII McEliece challenges including
  TII-252**, while itself stating the analysis "does not indicate HOVER would threaten Classic
  McEliece parameters in its present form".
- **Dated status line for the index:** as of 2026-08-30 the Classic McEliece **hold-out** dispute is
  still four-paper, with 1747 now at a 2026-08-30 revision adding a Frobenius-alignment method;
  separately 2026/1778 (HOVER) is a fifth McEliece cryptanalysis paper in the same wave on the
  higher-order-vanishing route. Both artifacts frozen in `mceliece/prim/` with the superseded v2
  copy **kept alongside rather than overwritten**.
- **UNIFORM $\delta$ RE-RUN LANDED, and the repair is a genuine unification.** Frozen at
  `campaigns/2026-08-30T14-09Z_57200ADD/`. All 13 instances through **one** code path, all exact:
  $\delta$ exact on 13/13, **zero aborts**, zero $L_i(a_i)$ failures across **6,144** unit checks
  ($64,64,64,128,128,256,256,512,512,512,1024,1024,2048$), zero off-diagonal violations in 2,600
  seeded probes, zero assembly spot-check failures in 39 rows; 0.02 s at $m=6$ to 285.91 s at
  $m=11$, 382 s total. The load-bearing degree condition $\max_j\deg f_j\le D$ was machine-checked
  per instance with **abort-not-warn** semantics, as required by the owner's counterexample ($f$ and
  $f+\Pi$ agree at every support point yet differ), and held on all 13.
- **The question that mattered: $m=6/7/8$ REPRODUCE under the new route.** Those rows had used the
  full $k\times n$ grid; they come out identical under the Lagrange reduction. So the repair
  **unifies** the ladder rather than substituting a weaker test that happened to agree — no finding,
  no escalation. **The index may now assert the uniformity it had been asserting; it is true as of
  this run.**
- **$\alpha$ was NOT relabelled, and the agent confirmed it had not seen the consequence.** The
  owner's amendment held: replacing the per-instance $\alpha$ measurement with a citation to Apon's
  Lemma 3 — the very lemma the measurement tests — would have moved the headline from *"Apon's
  mechanism confirmed by independent measurement"* to *"our data is consistent with Apon's lemma"*
  **without a single number changing.** $\alpha$ stays measured on all 13 rows; the artifact records
  that reasoning beside the label, and the $m=12$ policy (derive $\alpha$ from exact $\delta$ +
  Lemma 3 **only** if measurement is unaffordable there, labelled `CITED-DEPENDENCY` as a different
  *kind* of evidence) is recorded as a conditional **no row currently uses**.
- **Stopped rather than buying a 14th row.** $m=12$ not started, per the owner's instruction that a
  13-instance ladder uniform in $\delta$ with $\alpha$ measured throughout beats a 14th row bought
  by loosening what the ladder means.
- **Agent handoff, and it argues against the obvious next step.** The highest-value next unit of
  compute on this target is probably **not** $m=12$. It is measuring the **failure rate of Apon's
  $\Delta_{p,q}\ne0$ genericity hypothesis** over seeds and support orderings. The 13 rows observed
  $13/13$ nondegenerate — **but that is an incidental byproduct of instances chosen for the census,
  not an adversarial search**, so it is not evidence about the rate. The same in-sample/out-of-sample
  distinction that governs `rs-pe3d`'s pattern P applies here: a property that held on instances
  selected for other reasons has not been tested. **Neither disputant has that rate**, and Apon's
  §3.6 explicitly flags the all-$\Delta$-zero case as the hole in his theorem — so a *measured*
  degeneracy rate would be a genuinely new contribution to the four-paper dispute rather than another
  confirming row.
- **Two handoff facts so they are not reinvented:** the $m=12$ $\alpha$-guard fix is pre-registered
  in `pre_statement.md` ADDENDUM 2 but **unimplemented**, so a successor inherits the design *and*
  the label decision without the temptation to redesign either; and the family-block factorization
  ($\Psi=C\!\cdot\!J$, $\le9$ jet rows) in `census_small.py` with its regression counts is what
  makes any future scaling tractable.
- **Method promoted to a discipline rule (README 17a).** For sources without a dated API, track
  provenance by **content hash**, not by listing date: cache the artifact and compare size, hash and
  line count on every re-check. Strictly stronger than a date, because it detects a silent in-place
  revision that no date would reveal — as it did on the day it was adopted.

### Owner overclaim, caught and retracted the same day

- **`omega/` gate C's local-optimality conclusion is RETRACTED, and the error was the owner's.**
  The agent reported the weaker "certified-locally-optimal at $\pm10^{-7}$ within every feasible
  single-block family **tried**". The owner judged that this blurred the mechanism and pushed it
  to the universal form — "every single-block direction is infeasible", "the **unique** feasible
  point of the swept family", "local optimality by constraint rigidity … none exists" — which was
  then recorded in `RESULTS.md`. **That sharpening was wrong.** The agent implemented a bad
  instruction faithfully. Owner structural judgement now **0-for-6** this session.
- **The defect, precisely.** Six directions were tested — one $\pm$ repair pattern per region.
  Each `glob_dist` block carries **45** coordinates against **27** marginal rows plus **one**
  simplex row, so rank–nullity leaves a nullspace of dimension $\ge45-28=\mathbf{17}$.
  **Finitely many rejected directions cannot certify a positive-dimensional space.** And the
  supporting mechanism argument is circular in a subtle way: *"the marginals are linear in dist
  while dist$_{\max}$ is fixed, so any dist move breaks the marginal rows"* establishes only that
  a move **which changes the marginals** breaks them — a tautology. The step required is that
  **every** dist move changes them, which fails exactly when the marginal map has nontrivial
  kernel on the block, and the dimension count says it does. The repair pattern summed to zero by
  construction (hence cleared the simplex row bitwise) but was never built to lie in the
  marginals' kernel.
- **What survives:** the `rp` no-go, because there the dimension count genuinely closes — a
  *triple* pinned by three independent exact rows gives a **zero-dimensional** feasible set, with
  exactly-feasible $\lvert\Delta r\rvert\approx5.5\times10^{-17}$; the six probes' measured
  infeasibility (worst marginal move $7.000097472856237\times10^{-8}$ against a
  $1.1\times10^{-9}$ gate, $63.64\times$, owner-recomputed); the absence of any box below
  published $2.37155181$, so no record claim ever existed; Part 1's precision-invariance result,
  which depends on none of this; and both accessor fixes including bug 2's three-instrument
  closure at $3.5258684860650646\times10^{-11}$.
- **Status corrected:** gate C **PARTIAL**, not COMPLETE. Scoped claim of record — no feasible
  improvement among the six pre-registered single-block probes at $\pm10^{-7}$, and region_prop
  admits no feasible motion at all. Whether `glob_dist` admits a feasible improving direction is
  **OPEN**; settling it requires computing the equality/Jacobian nullspace on the block and
  certifying every remaining direction, **not** testing more hand-chosen patterns. Sampling more
  directions would repeat the error being corrected.
- **The lesson, and it is about the owner.** The same session caught three shared-instrument audit
  failures in agents' work — `kg/`'s finite-difference check against the analytic derivative of
  its own wrong function, `mm3/`'s three solvers over one modelling premise, `omega/`'s accessor
  comparing every row to a single $b$ — and then committed a dimension-counting error of the same
  family: accepting "six directions rejected" as evidence about a **space**, because the mechanism
  story sounded structural. **A mechanism that sounds structural is not a proof that the mechanism
  is exhaustive.**
- **Exact rank supersedes the owner's estimate, and the estimate's flaw was the assumption the
  owner had told the agent to check.** Agent SVD count: 45 coordinates; the 27 margin rows are
  three full-rank 9-row blocks whose spans **overlap**, so stacked margin rank is **24**; and the
  simplex row of ones lies **inside** the margin span, adding zero rank. Independent rows **24**,
  nullspace **exactly 21**. The owner's $45-28=17$ assumed independence — the very assumption the
  owner instructed the agent to verify and did not verify itself.
- **A symptom of rank deficiency was read as a symptom of rigidity.** The repaired probes cleared
  the simplex row *bitwise*, which both owner and agent took as the repair working. It was a
  **trace of the dependency**: the simplex row lies in the margin span, so anything satisfying the
  margins satisfies it for free. The kernel is simplex-orthogonal to machine precision
  ($\max\lvert v\cdot\mathbf1\rvert=5.55\times10^{-16}$).
- **The OPEN question is genuinely open.** All 21 kernel directions are two-sided at $p^*$
  (support min $1.28\times10^{-7}$ exceeds the margins), so the feasible set inside the box around
  dist is $\ge21$-dimensional; and the objective is **not** constant there — marginals invariant,
  but $H(\mathrm{dist})$, $p_{\rm comp}$ complete-split terms and the `part_frac` flow do vary,
  with float sensitivity at $\epsilon=10^{-7}$ spanning $-2.5\times10^{-9}$ to
  $+6.8\times10^{-9}$, **both signs**.
- **Owner-computed scoping that bounds the whole campaign in advance:** at that sensitivity a
  fully successful certification improves our enclosure by $\approx10^{-9}$ **at most**. Against
  the $2.02585446\times10^{-6}$ gap to published, the shortfall factor is $\mathbf{810.34}$ and
  reaching published would need $\epsilon\approx8.103\times10^{-5}$ — three orders outside the
  $10^{-7}$ box. **The value is the certification, not the number**, and the pre-statement must
  say so in its first paragraph so no reader mistakes the outcome for a record attempt.
- **`omega/sms` checked rather than assumed — second branch as well.** No linear equality row
  touches $\omega$ or `single_mat_size` ($K$ pinned by bounds, outside box scope); the only acting
  constraint is the Schönhage line $\mathrm{sms}\cdot\omega=T-R$, **one** equation on **two** box
  coordinates, codimension 1, so **not frozen** — but total feasible $\omega$ excursion along it
  is $\approx4.5\times10^{-8}$, under half the box radius, and real motion needs the compensating
  $R$ already marked OPEN. No free second no-go, no live single-coordinate descent.
- **Rule adopted: OWNER-SUPPLIED FIGURES ARE INPUTS, NOT EVIDENCE.** Twice today an owner number
  moved through an agent unchecked — the König "cap" hypothesis (pushed back on only after a
  full-text read) and the $27/45/1$ row counts (corrected only when a rank computation was
  mandated). Any claim resting on an owner figure must be independently re-derived before it is
  written. **The rank/dimension check is now the required instrument for every "frozen by rows"
  claim in this repository.**
- **Full 21-dimensional subspace certification APPROVED as the real gate C**, under four
  conditions: the $810\times$ shortfall scoping stated in the pre-statement's opening paragraph;
  coverage of the whole subspace or an explicit failure-to-certify (no shrinking to directions
  where the remainder happens to bound); a **hard one-hour checkpoint** on the Hessian-Lipschitz
  remainder, after which the scoped negative stands as final; and escalation before any wording
  on either outcome. Route: interval directional derivatives along the 21 exact SVD directions
  plus a certified box-wide Hessian-Lipschitz remainder (entropy Hessian exactly $-1/t$ per
  coordinate, $p_{\rm comp}$ low-degree).
- **Certification BLOCKED before compute: the kernel basis was a float SVD basis.** Two
  consequences, the second fatal. (i) `kernel_basis_V.npy` is only *approximately* in the kernel —
  simplex-orthogonality holds to $5.55\times10^{-16}$, i.e. float noise, not zero — so
  $p^*+Vc$ satisfies the 27 marginal rows to $\sim10^{-16}$ rather than exactly. That is inside
  the $1.1\times10^{-9}$ program tolerance and perturbs the objective by $\sim10^{-25}$ against a
  $1.06\times10^{-8}$ linear swing, so it is a **labelling** problem: nothing built on it can carry
  `MACHINE-VERIFIED`. (ii) **Fatal: "three zero singular values" does not prove exact rank 24.**
  The risk is asymmetric — if the true exact rank is *higher*, $\mathrm{span}(V)$ merely contains
  infeasible directions and a cap only gets stronger; but if it is *lower*, the true nullspace is
  **larger than 21** and $\mathrm{span}(V)$ **misses a kernel direction**, making a cap over
  $\mathrm{span}(V)$ silent about it. **So the count that superseded the owner's estimate is
  itself unproven** — the owner's $17$ assumed independence it had not checked; the agent's $24$
  rests on a float rank decision. Same class of error.
- **Fix ordered, and it is cheap:** extract the exact $27\times45$ coefficient entries
  (rationalize exactly and document any entry that is not exactly representable), compute exact
  rank and nullspace over $\mathbb Q$ by fraction-free elimination, compare against the SVD's 24
  and correct the pre-statement per rule 5 if they differ, verify $A\cdot V_{\rm exact}=0$
  **identically** and — separately — verify **exactly** that the simplex row lies in the margin
  span, since that claim also came from the float rank computation. Interval certification is an
  acceptable alternative provided enclosures exclude zero for the nonzero singular values and
  contain zero for $A\cdot V$; a tolerance comparison is not certification. The 60-minute remainder
  checkpoint restarts only once the basis is exact.
- **A second scoping trap closed in advance rather than after.** The feasible set in coefficient
  space is a **polytope**, not a cube: the coordinate box $\pm10^{-7}$ pulls back through $V$ to a
  polytope whose inscribed cube has side $r_c=3.3631\times10^{-8}=10^{-7}/2.97$, the $2.97$ being
  the row-sum of $\lvert V\rvert$. Along directions with small entries a single coefficient may
  range far beyond $r_c$ while staying in the box, so
  $C\subsetneq(\text{feasible kernel region}\cap\text{box})$ strictly. **A cap over $C$ is a cap
  over $C$**, and "certified: no improvement over the 21-dimensional feasible subspace at
  $\pm10^{-7}$" would be false in exactly the way the retracted claim was false. The result
  sentence must name its certified domain and, if $C$, state the uncovered remainder explicitly;
  "certified" without a named domain is how the previous overclaim happened. Re-centering on a
  certified linear minimizer is permitted but may **never** convert a failure-to-certify into a
  success — moving the domain after seeing the outcome is domain-shopping, the same failure mode as
  shifting a search box to find a win, and any re-centered result is a separate statement with its
  own named domain reported alongside the original attempt.
- **Third instance today of one shape, now the session's central lesson.** An approximation stood
  in for the set it represents, three times: **six probes** for a 21-dimensional space, an
  **inscribed cube** for its circumscribing polytope, a **float SVD basis** for an exact kernel.
  **The instrument must be able to establish coverage, not merely be consistent with it.**
- **The rank-24 deficiency is now EXPLAINED structurally, and the explanation cross-checks the
  rref by an independent route.** Ordered the agent to name the three left-nullspace relations,
  because an unexplained rank deficiency could equally be a duplicated or mis-transcribed row —
  which would have been a third transcription bug. All three are **structural**, no duplicates
  (27 rows pairwise distinct, checked): **dep0**, the shared grand total
  $-\sum_j f_{0,j}+\sum_j f_{1,j}=0$; **dep1**, a first-moment relation
  $\sum_j-(8{+}j)f_{0,j}+\sum_j(8{-}j)f_{1,j}+\sum_j(8{-}j)f_{2,j}=0$ with exact size check
  $\sum(8{+}j)(9{-}j)=480=\sum(8{-}j)\cdot2(9{-}j)$; **dep2**, its variant, $435=435$ exactly.
  Owner re-verified both size checks.
- **Owner catch: the agent's accounting listed FOUR relations where rank 24 permits exactly
  THREE.** dep0, dep1, dep2 were named, and the second total relation
  $t_0\text{-total}-t_2\text{-total}=0$ was separately verified — four in total, against a
  left-nullspace of dimension $27-24=3$. Resolved symbolically by the owner: adding the moment
  relations coefficient-wise in $j$ gives $f_0\mapsto-1$, $f_1\mapsto0$, $f_2\mapsto+1$, i.e.
  $$\text{dep1}+\text{dep2}=-\textstyle\sum_j f_{0,j}+\sum_j f_{2,j}=0,$$ **exactly the second
  total relation.** So $\{$dep0, dep1, dep2$\}$ is a basis and the fourth relation is redundant.
  The count is therefore consistent with the proven rank by an argument **independent of the
  rref** — a second instrument agreeing with the first.
- **This also retires the owner's generic-model discrepancy properly.** The owner's toy
  three-block marginal model predicted rank 25 (two total relations, nullspace 20) against the
  actual 24. The resolution is not "one extra dependency" but a **different decomposition**: it is
  *one* total relation plus a **2-dimensional moment family**, and the model's second total
  relation lies inside that moment span. **The laser-method margins are moment-coupled, not merely
  total-coupled** — a program-structure fact worth more than the count that prompted it.
- **Certification re-routed on an owner observation: the inscribed cube was a basis artifact.** The
  exact rational rref basis reported $r_c=4.7619\times10^{-9}$ with "margin factor $21.0$" — and
  $10^{-7}/21=4.7619\times10^{-9}$ *exactly*, so the margin factor **is** the nullspace dimension.
  That is the signature of an **unnormalized** rref basis ($\pm1$ pivot entries, row-sums of
  $\lvert V\rvert$ near the dimension), not a property of the subspace. Certifying that cube would
  have covered $1/7.06$ of the already-scoped SVD-basis domain, with the linear swing falling from
  $\pm1.06\times10^{-8}$ to $\pm1.50\times10^{-9}$, best-case improvement to
  $\approx1.90\times10^{-10}$, and the shortfall factor to published rising from $810$ to
  $\mathbf{10\,677}$ — while leaving nearly the whole feasible region uncertified.
- **Better route adopted: certify the polytope, not a cube.** Parameterize by the coordinate
  perturbation and take the domain to be the natural set
  $D=\{\delta: A\delta=0,\ \lVert\delta\rVert_\infty\le10^{-7}\}$ — exactly kernel $\cap$ box, with
  no inscribed-cube loss and no basis dependence. The linear term is then an **LP** (45 variables,
  24 equality rows, box bounds) bounded rigorously by exhibiting a dual-feasible $(y,\mu^\pm)$ and
  verifying weak duality in interval arithmetic — no sampling anywhere. The Hessian remainder over
  the full coordinate box is the **same** computation as over a sub-cube, since the kernel
  constraint can only shrink the achievers, so the stronger claim costs nothing extra. Caution
  issued: the dual point must be rounded to exact rationals and verified **dual-feasible** exactly,
  since a merely approximately-feasible $y$ gives a non-rigorous bound — a slightly suboptimal but
  exactly feasible dual still yields a valid bound, and **a valid loose bound beats an invalid
  tight one.**
- **Agent's lesson line, adopted repo-wide:** *a round number in a scaling factor is a basis
  artifact until proven geometry.* It generalizes past bases — a suspiciously round constant
  anywhere in a pipeline is a property of the representation until someone proves otherwise.
- **Certification attempted and it FAILED TO CERTIFY — reported as such, not patched, not
  shrunk.** The naive two-sided certified enclosure over the enclosing box gives
  $\omega\in[2.3715431556360005,\ 2.3715624359361276]$, width
  $1.9280300127\times10^{-5}$, against a first-order LP signal over $D$ of
  $\min G\cdot\delta=-1.5816497000997742\times10^{-7}$ — the enclosure is **$121.9\times$ looser
  than the signal** (owner-verified). Subdivision cannot fix it: widths are linear in radius, so
  closing the gap needs radius $\sim10^{-10}$, i.e. $(10^3)^{45}=10^{135}$ boxes. What is missing
  is a **certified gradient/slope enclosure** through the transcribed tree plus a second-order
  remainder — new machinery, correctly outside the 60-minute budget. The box was verifiably
  **live** this time (point-eval width exactly $0.0$; glob `num_block` widths
  $6.656/6.616/6.616\times10^{-6}$ for region 0 against $3.614\times10^{-20}$ untouched).
- **The polytope route also raised the POTENTIAL by two orders of magnitude, not just the rigour.**
  Best conceivable gain over the true domain $D$ is $1.58\times10^{-7}$, i.e. **$12.808\times$**
  short of the $2.0258544615181506\times10^{-6}$ gap to published — against $810\times$ for the
  SVD-basis cube and $10\,677\times$ for the rref cube. The question moved from *four orders
  short* to *one order short*: from a dead end to a live open problem.
- **The agent's Lagrangian instrument is strictly stronger than the caution the owner gave.** The
  owner warned that an only-approximately dual-feasible $y$ yields a non-rigorous bound. The
  concern evaporates: since $A\delta=0$ identically on the kernel, $y^{\!\top}A\delta=0$ for
  **any** $y$, so $G\cdot\delta=(G-A^{\!\top}y)\cdot\delta$ and Hölder gives
  $$G\cdot\delta\ \ge\ -R\lVert G-A^{\!\top}y\rVert_1\quad\text{for arbitrary }y,$$
  with no dual feasibility to verify — optimizing $y$ merely tightens it. With the LP's $y$ the
  bound matches the LP optimum to $1.6\times10^{-22}$; on a $10^{-6}$ rational grid it costs
  $1.25\times10^{-12}$ and stays rigorous. Recorded as the reusable instrument for whoever builds
  the gradient enclosure.
- **FOURTH instance of the session's failure shape, and the most transferable: a solver whose
  tolerance equalled the quantity being checked.** `scipy linprog(method="highs")` with variable
  bounds of magnitude $10^{-7}$ returned a point **violating those bounds by 100%** (two
  coordinates at $-2\times10^{-7}$ against $-10^{-7}$) and an objective **14.12% below** the true
  optimum — because HiGHS's **default primal feasibility tolerance is $10^{-7}$**, exactly the box
  radius. The solver was not wrong by its own lights; it answered a different question, plausibly.
  Fix: **non-dimensionalize** ($\delta=Ru$, $\lVert u\rVert_\infty\le1$) rather than tightening
  tolerances and hoping — bound violation then exactly $0.0$, kernel residual exactly $0.0$. The
  agent caught it only because it re-verified the returned point against its own constraints;
  earlier unscaled values ($-5.6\times10^{-8}$, $-1.805\times10^{-7}$) are discarded per rule 5 and
  were never recorded. **Broadcast to `Mm3GateC`** (its gate B used HiGHS ILP `INFEASIBLE` as one
  of three routes, and an `INFEASIBLE` that is really *infeasible within tolerance* is the mirror
  image of the Heule–Kauers–Seidl Challenge-2 defect) **and to `DelcapGateC`** (BA locating path
  and finite-length LP work, resolving $10^{-11}$–$3\times10^{-8}$). **Repo rule now: re-verify
  every solver's returned point against its own constraints in independent arithmetic — a solver's
  reported objective and its reported feasibility are both claims, not evidence.**
- **Gate C closed. The AD build is NOT approved**: new multi-hour machinery whose total success
  caps out $12.8\times$ short of published, so it cannot yield a record and would tighten only our
  own enclosure. Now precisely specified for a future campaign, which is the valuable part. A clean
  complete negative beats an extended incomplete positive.

### `rs-pe3d/` gate B in flight — beta is EXCLUDED as the driver, and a discriminating test was pre-registered

Pre-statement committed **before** any computation at
`campaigns/2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB/pre_statement.md` (18.5 KB, verbatim
Def 4.1 / Conj 4.2 / Thm 2.1 / Thm 4.11 / Lemma 5.4 balance quotes from the agent's own PDF read,
the $\beta$ definition with proof, the 13-row grid, budget, adjudication rules, rule-7 statement).
$\beta(s):=\max_i s_i/\min_i s_i$, proved to be the minimal $\beta$ satisfying the conjecture's
hypothesis and to specialize to Thm 4.11's $K$ at $r=2$. **12 of 13 rows exact and complete**; only
B11 $=(241,(3,5,16),(1,1,1))$, $\beta=16/3$, in flight (2,304,200 supports, 45.6% done after a
$3.7\times$ vectorized mod-$q$ rref rewrite). Anchors reproduced with the sweep code as run: the
unit $(5,(2,2,2),(1,1,1))\to\rho=1/3$ exact over $5^7-1=78{,}124$ nonzero vectors, **plus a
structurally independent second route** (V-basis coefficient enumeration + closure) also giving
$1/3$; and all seven frozen gate-A rows matched exactly, zero mismatches.

- **$\beta$-as-driver is not merely unsupported — it is EXCLUDED by the completed rows.** At $q=31$
  the ratios are **interleaved** in $\beta$: $\beta=5/2\to3/5$, $\beta=3\to1$, $\beta=5\to3/5$. A
  value of $1$ sits strictly between two values of $3/5$ as $\beta$ increases, which no monotone
  function of $\beta$ can produce. Same shape at $q=241$: $\beta=5/3\to1$, $\beta=5/2\to3/5$,
  $\beta=8/3\to1$. This is the contract's **third** admissible outcome ("neither"), which was
  explicitly permitted rather than forced into the graceful/collapse dichotomy.
- **A sharper pattern, owner-tabulated and holding on all 12 rows:**
  $$\{2,3\}\subseteq\operatorname{orders}(s)\iff\rho=3/5 .$$
  The five $3/5$ rows all contain both 2 and 3; the seven non-$3/5$ rows ($\rho=1$ four times,
  $1/2$, $3/8$, $1/4$) contain neither pair. So the driver appears to be the triple's **arithmetic
  shape**, not its imbalance.
- **The in-flight row DISCRIMINATES, and the prediction was pre-registered before it landed.** B11's
  orders $(3,5,16)$ contain 3 but **not** 2, and its $\beta=16/3$ is the largest in the grid. So the
  arithmetic-shape pattern predicts $\rho=1$, while any $\beta$-degradation story predicts the
  **lowest** ratio in the grid — opposite predictions on one row. The agent was directed to commit
  both predictions as a dated addendum **before** the census completes: *a pattern noticed after the
  datum arrives is a story; the same pattern committed before it arrives is a test.* Owner
  calibration attached to the suggestion in writing: hypotheses stand at **0-for-8**, so the
  $\{2,3\}$ pattern is to be killed rather than confirmed, and a third value (neither $1$ nor $3/5$)
  would falsify **both** and be the most informative outcome available.
- **The most publishable finding in the target, and the agent had buried it.** The **AI candidate
  proof itself carries K-comparability** — $\rho=(1/(2dK^d))(\varepsilon/(C_0dK)^{10d})^{3d\,6^{d-2}}$
  — as do Thm 4.11 ($K^{-1}\le n_1/n_2\le K$) and Lemma 5.4 (4-balancedness). **Only the informal
  Conjecture 1.4 drops the balance hypothesis.** That is exactly the outcome gate B's contract named
  as its best case: *a precise, publishable correction to the conjecture's statement rather than to
  its truth*. The gap is not between a proved theorem and a bolder conjecture, but between the
  formal statements — which all carry balance — and one informal restatement that omits it. Ordered
  into its own section with the three quotes side by side and line numbers; it stands independently
  of B11.
- **B11 LANDED AND THE PRE-REGISTERED TEST RESOLVED: branch 1.** $\rho(\mathrm{B11})=1$ **exactly**
  at $(241,(3,5,16),(1,1,1))$, $\beta=16/3$ — witness weight 3, $\delta=3$, closure exact with all
  cheaper line-tuples absent, independently re-verified in exact $\mathbb F_q$. Census complete: all
  2,304,200 weight-$\le3$ supports, 80 candidates, 3448 s. **13 of 13 rows landed**, all exact
  rationals, 12 re-run in a fresh process with zero ratio disagreements.
- **$\beta$-degradation is FALSIFIED AT THE OPPOSITE EXTREME, not merely unsupported.** It predicted
  $\rho(\mathrm{B11})<1/4$ — strictly *below* the grid minimum, since $\beta=16/3$ is the most
  imbalanced point swept. The observed value is $1$, the grid **maximum**. Predicted floor, observed
  ceiling. Owner-verified: at $q=241$ in increasing $\beta$, $5/3\to1$, $5/2\to3/5$, $16/3\to1$ —
  $\beta$ grows $3.2\times$ with $\rho$ returning to $1$ while the *intermediate* $\beta=5/2$ sits at
  $3/5$. No monotone function of $\beta$ does that.
- **Verdict: the third outcome, NEITHER graceful degradation nor collapse at a $\beta^*$.** Exclusion
  **E is a conclusion**, independent of B11 and strengthened by it: on the swept set
  $\rho^{\rm window}$ **is not a function of $\beta$ at all**. There is no $\beta^*$ to locate, so no
  threshold claim and nothing to escalate as a hypothesis correction; the grid-resolution statement
  applies vacuously.
- **Pattern P is EVIDENCE, not a theorem, and the in-sample/out-of-sample split is stated.** All 13
  rows partition perfectly — 5 contain both 2 and 3 and all five are $3/5$; 8 do not and none is
  $3/5$. In-sample separation under a random-labelling null is $1/\binom{13}{5}=1/1287=7.8\times
  10^{-4}$, **and that figure is in-sample**: 12 of the 13 rows are the ones the pattern was found by
  inspecting. What P has is a perfect fit on 12 in-sample rows plus **one out-of-sample
  pre-registered confirmation** carrying a cryptographic order proof. Recorded as *"P survived its
  first prospective test"*, never as established, with its falsifier named (any triple containing
  both 2 and 3 with $\rho\ne3/5$, or any lacking the pair with $\rho=3/5$) and its limits stated —
  P is a statement about **orders** at $t=(1,1,1)$, and the three $q=13$ rows show $\rho$ varying
  with $t$ at fixed orders ($1/2$, $3/8$, $1/4$), so P cannot be the whole story.
- **Two further owner checks on the frozen 13-row table, which carries `wt`, `delta` and
  `pairwise_coprime` — variables not available when P was proposed.** (i) $\rho=\mathrm{wt}/\delta$
  in every row, confirming the definition is applied consistently. (ii) **P is not a shadow of
  coprimality:** it holds *separately within both classes* — among the 7 pairwise-coprime rows the
  three $(2,3,5)$ instances give $3/5$ while $(3,4,5)\times2$, $(3,5,8)$, $(3,5,16)$ give $1$; among
  the 6 non-coprime rows $(2,3,10)$ and $(3,10,2)$ give $3/5$ while $(2,2,4)\times3$ and $(2,6,5)$
  do not. Since coprimality is a hypothesis TR26-150 itself carries (line 861, *"pairwise coprime,
  4-balanced orders"*), ruling it out as the driver matters.
- **A separate observation, labelled as such and not a claim:** on these 13 points **coprimality
  appears to restrict the ratio set**. The 7 pairwise-coprime rows take only **two** distinct values,
  $\{1,3/5\}$; the 6 non-coprime rows take **five**, $\{1/4,3/8,1/2,3/5,1\}$. That is an observation
  on 13 in-sample points with no prospective test behind it — recorded as a candidate for the next
  campaign's pre-registration, explicitly **not** as a finding, and it is the natural companion
  hypothesis to P since both concern arithmetic structure rather than balance.
- **Owner calibration:** P was the owner's suggestion, tabulated from the agent's data and handed over
  explicitly flagged as a hypothesis to kill rather than confirm, with the owner's score at 0-for-8
  at that moment. It survived one test — the owner's first structural suggestion of the session that a
  measurement did not immediately destroy, and one survival is not a record, which is precisely why
  the wording above is hedged.
- **OWNER VERIFICATION FAILURE, the fifth instance of rule 14 and the first committed by the owner
  while enforcing it.** The AI candidate's constant is
  $\rho=(2dK^d)^{-1}(\varepsilon/(C_0dK)^{10d})^{d/\gamma_d}$. The agent first reported the $d=3$
  degradation as $K^{-4863}$; the owner "confirmed" it by inferring that the ASCII exponent `3d`
  must mean the tower $3^d$ **because that reading reproduced the agent's number**. Both were wrong.
  The proof pins the exponent by its own definitions — line 245 sets $\gamma_1=1$,
  $\gamma_e=1/(3\cdot6^{e-2})$; line 489 gives $\gamma_e=\gamma_{e-1}/6$; line 933 states the
  exponent *is* $d/\gamma_d$ — so in exact `Fraction` arithmetic $d/\gamma_d=6,54,432,3240$ at
  $d=2,3,4,5$, matching the **linear** $3d\cdot6^{d-2}$ at every $d$ and excluding the tower
  ($9,162,2916,52488$). Hence exponent $54$ at $d=3$ and $\rho\sim K^{-1623}$, **not** $K^{-4863}$.
  **The owner's check was circular:** it inferred the reading from a quantity derived out of the very
  formula in question, so it was guaranteed to reproduce whatever reading generated the number being
  checked — two apparently independent confirmations of one wrong figure, one derived from the other.
  The agent caught it by re-deriving from the definitions, retracted its own figure inline with its
  cause, **and recorded that the owner's confirmation had matched the error rather than catching it.**
  Corrected in `RESULTS.md` with both readings' values at $d=2..5$. **Owner checks now 0-for-9.**
  Rule adopted for this target: when a formula's reading is in doubt, resolve it from the source's
  **definitions**, never from a number computed with it.
- **And a correction to this repository's own framing, forced by the agent's first-hand read.**
  `RESULTS.md` had asserted since gate A that the proved 2-D Theorem 2.1/4.11 carries a balance
  hypothesis *"that Conjecture 4.2 drops for higher $r$"*. **Conjecture 4.2 does not drop it** —
  lines 740–748 carry *"balance parameter $\beta\ge1$ … $\beta^{-1}\le s_i/s_j\le\beta$"*.
  Balance is present in Conjecture 4.2, Theorem 4.11 (1280–1285), Lemma 5.4 item 2 (1632), Section 6's
  own summary of its hypotheses (1604), **and in the unverified AI proof** (Theorem line 8, Prop 5
  line 252); it is absent **only** from Conjecture 1.4's informal phrasing (224–226). Retracted
  inline per rule 5 with the error attributed to us, not the paper. **This also undercuts gate B's
  original premise** — the gate was built on the conjecture omitting a hypothesis the theorem has,
  which is false for Conjecture 4.2 — so the defensible result is the provenance finding itself: the
  correction belongs to one informal restatement, and Conjecture 4.2 needs none on this axis.

### `oct-rank/` gate C COMPLETE — the last "nothing filed" item in the repo is closed

Campaign `2026-08-30T12:15:25Z_gateC4_c598b05f…_0481d816063f` (agent `OctRankGateC`, hashes
verified). All three subtasks filed; **no refutation, no discrepancy with the paper, no escalation
triggered.** One in-flight witness-orientation bug was caught by the exact checks, fixed, and
recorded as a rule-5 retraction in `VERDICT.md`.

- **$\mathrm R_{\mathbb R}(\tau)=7$, `MACHINE-VERIFIED`.** Upper: Krawczyk replay re-run PASS,
  exact $K<1$, margins $9.946\times10^{-6}$. Lower: every consumed polynomial identity re-derived
  exactly (sympy $\mathbb Z[\dots]$ + `fmpq` grid, **no floats**), on a `CITED-DEPENDENCY` Lean
  gate-A skeleton.
- **Sharpness of $\tfrac52n-2$: $n=2$ is now FULLY SELF-CONTAINED.**
  $\mathrm R_{\mathbb R}(T_{\mathbb C})=3$ `MACHINE-VERIFIED` in **both** directions — an exact
  three-term witness, plus the target's own exact lower proof: det-multiplicativity forces
  $x^2+y^2=W\cdot(\text{two linears})$, and coefficient matching yields
  $(W f_{10} f_{21})^2=-1$ over $\mathbb R$, with an exact Sturm zero-root certificate for
  $t^2+1$. **Owner-verified independently:** $x^2+y^2$ is irreducible over $\mathbb R$ and factors
  only as $(x-iy)(x+iy)$ over $\mathbb C$; $t^2+1$ has zero real roots (discriminant $-4$); and
  $z^2=-1$ has empty solution set over $\mathbb R$. The argument is sound, and $3$ agrees with
  Gauss/Karatsuba above and with `2608.16649`'s stated rank $3$ for $\mathbb C$, $8$ for
  $\mathbb H$. For $n=4$: $\mathrm R_{\mathbb R}(T_{\mathbb H})=8$, exact eight-term witness above,
  lower via peel(2)+pencil(6) as `CITED-DEPENDENCY`.
- **Gate C4 — the $n=8$ peeling-constant probe: NEGATIVE, as pre-registered.** This item had stood
  as `HUMAN-AUDIT-PENDING` scratch with **nothing filed**; it is now closed. All probes exact
  (`fmpq`/sympy, no floats): pencil rank $(1,C)=12$ **exactly** on the consumed class, via an
  entrywise-verified 12-term witness for $(I_8,J_8)$ direct-summed from the verified $2\times2$
  Gauss witness, plus the Lean floor 12; **every** stopping point $j=0..6$ gives exactly $18$; and
  the `key_bound` inequalities are exactly tight ($UD=0$, $UA=I$, $UB=J$, $\operatorname{rank}D=4=n/2$,
  $\dim\ker U=4=r-n$). **So $18$ is saturated inside the substitution/peeling + Thm-2 machinery and
  no parametric choice lifts it.**
- **Rule-7 boundary, stated:** not swept and therefore open — rank $\ge14$ for 3-slice
  $L$-families $(L_u,L_v,L_w)$, which *would* lift $18$ to $19$; and **all non-peeling routes**.
  The negative is about the peeling machinery, not about $\mathrm R_{\mathbb R}(T_{\mathbb O})>18$.

### `kg/` gate C COMPLETE — a certified cap on the septic family, and a ceiling that reframes the $b_5$ route

Frozen at `campaigns/20260830T130130Z_05F93E5B_83eb3752116e/` (REPORT.md, pre-statement with Addenda
1–4, code, logs, manifest, 31-line `sha256.txt` including the external `CITED` inputs). README
appended below the contract; ladder untouched. **No paper-side disagreement anywhere** — every
shared quantity agrees ($b_1=0.8815738220495995485818877$ vs their $\ge0.881573822049$; head
$1.132885992768967\times10^{-5}$ vs their $\le1.1328860\times10^{-5}$; $11/6$; $6\pi/11$).

- **Anchor and tripwires first, reported as numbers.** avg$(S)=123.37760876657072$ (rel
  $1.9\times10^{-9}$ against the owner's frozen $123.377609$), $B_3=14.245983003864794$ (rel
  $2.7\times10^{-10}$). T1 $\pi$-periodicity $9.93\times10^{-16}$; **T2 reflection worst
  $7.94\times10^{-15}$ against the frozen $7.936\times10^{-15}$ — no degradation** across 24
  node-pair equalities; T3 $\min S=1.266545$ exactly at $\pi/2$, $\max S=310.80695994744156$. The
  right generalization was made: exponents $1,3,5,7$ are all odd, so $\rho_j(-t)=-\rho_j(t)$ is an
  **identity** for the whole septic family and both tripwires extend without re-derivation.
- **(a) CERTIFIED CAP: the septic family does not beat $3.47\times10^{-4}$.** The best certified
  point in the entire sweep is the **paper's own** $(s_3,s_5)$ at $s_7=0$. Pruning bound
  $\gamma^*\le b_1(V_{\min})-\text{head}_{\rm lo}$ needs no tail (tail $\ge0$, $\lvert b_m\rvert\ge0$).
  Depth-7 pre-registered: $99.97811\%$ pruned, 2409 evaluations; depth-10: $99.998739\%$, survivor
  hull $s_7\le0.019922$; anisotropic refinement of the provably-containing survivor slab:
  $99.9999\%$ pruned, residual $3.088\times10^{-11}=1.43\times10^{-8}$ of the box. Inside it the
  septic-natural triple annihilation $b_3=b_5=b_7=0$ at $s_7=0.00036449$ **loses** by
  $-3.37826580833\times10^{-8}\pm3.99\times10^{-20}$.
- **The mechanism is the valuable part, and it is certified.** $\mathrm db_1/\mathrm dV=-0.787793605140965$;
  $b_7=1.17615565426\times10^{-7}$ dies at $u_7=1.32852960100\times10^{-7}$, giving a
  $\lvert b_7\rvert$ gain rate of $0.885306321646$ per unit $u_7$. **The gain rate beats the cost
  rate but saturates at the zero-crossing while $m\ge9$ keeps growing**, so only
  $5.36\times10^{-8}$ of $1.18\times10^{-7}$ survives against a linear $8.74\times10^{-8}$ cost.
  That explains why no refinement of the search can help — strictly better than "the search found
  nothing". Ten monotone certified continuation nodes on the paper's own $b_3=b_5=0$ manifold
  ($s_7=0.0002\to-8.43\times10^{-9}$ … $s_7=0.05\to-5.04\times10^{-3}$) and Nelder–Mead converging
  to $s_7=0$ from both starts close it.
- **Addendum-3 catch: a convention defect that would have silently voided every dependent run.**
  $A_{a,b}$ depends on $\vartheta=\eta/\sqrt V$ (paper lines 1717–1719, 685), so an $\eta$-fixed
  septic sweep moves $\vartheta$ as $V$ changes and **invalidates the imported A-grid**. Found
  *before* the dependent runs; fixed $\vartheta=\vartheta_{\rm paper}=0.128957369921412$ via
  $\eta(s)=\vartheta_{\rm paper}\sqrt{V(s)}$ ($\eta$ is free per line 675); anchor unaffected at
  $s_7=0$. **Cost stated as scope, not buried: the $\vartheta$ axis is NOT swept** — it needs
  certified A-grid regeneration, 16,002 one-dimensional integrals, priced $\le8$ h, not run.
- **A better headline declined for the right reason.** The corrected float $B_3$ gives
  $K_G\le1.781841202333086968398$ against the cited-tail chain's $1.781841328136936456209$ —
  **stronger by $1.25804\times10^{-7}$, owner-verified** — but float quadrature sits in that chain,
  so it is `COMPUTATIONAL-EVIDENCE` and was **not** promoted. Label discipline operating against the
  repository's own interest, which is the only circumstance in which it means anything.
- **(b) The degree-5 affine program: a reduction, a certified ceiling, and a provable no-go.** With
  $A_\lambda=2b_1-b_3-\lambda b_5\le C_\lambda$ (affine, mixture-stable) and
  $a_5=(3b_3^2-b_1b_5)/b_1^7$, the barrier becomes $\gamma\le\min(11/12,\,C_\lambda/2)$, so the
  $b_5$ route helps **iff** $C_\lambda<11/6$, giving $K_G\ge\pi/C_\lambda$; $\lambda=0$ reproduces
  the paper exactly. Owner-verified certified values: the hyperplane pair has $b_1=1$, $b_3=1/6$,
  $b_5=3/40$ and $A_0=2-1/6=11/6$ **exactly, so $11/6$ is sharp**; the Naor–Regev ceiling forces
  $C_\lambda\ge\pi/K_G^{\rm upper}=1.76311583499934371781$, whence
  $\lambda^*=0.9362333111198615$.
  * **A genuine no-go:** at $\lambda=1$, $C_1=11/6-3/40=1.7583333333333$ forces
    $K_G\ge1.786687765074764$, which **exceeds our upper bound** — so a hyperplane-attained
    $b_5$-strip is **impossible** for every $\lambda>\lambda^*$.
  * **The admissible range is $0<\lambda\le0.9362333111198615$**, and certifying $C_\lambda<11/6$
    anywhere in it improves the lower bound.
  * **At the top of that range the improvement would close the ENTIRE remaining gap:** at
    $\lambda^*$ the barrier gives $\gamma\le0.881557917499671859$ and
    $K_G\ge1.781841323880437288$. **Owner precision note:** that this equals our upper bound to
    $\pm3.59\times10^{-90}$ is **true by construction** — $\lambda^*$ is *defined* as the $\lambda$
    at which the barrier's implied lower bound meets our upper bound — and must not be read as an
    independent coincidence. The content is that **the degree-5 program's admissible range extends
    exactly to the point where it would determine $K_G$**, so it is not a small refinement of the
    degree-1-and-3 program; its endpoint is equivalent to the full problem.
- **The agent proved its own swept set is the wrong set.** Families F1 (1-D strip, 11 members) and F3
  (2-D Davie–Reeds, 13 members) in closed-form Arb: $b_5$ changes sign inside both; the best $A_0$
  with $b_5\le0$ is $1.47033139820178559$ ($0.363$ below $11/6$); and the swept
  $\lambda_{\max}=9.14497596448671$ **far exceeds** $\lambda^*$, which **proves the swept families
  do not contain the degree-5 extremisers** — independently matching the paper's own line-1013
  remark that stronger strips need $\lVert P_1k\rVert^2$ control. A proof that one's own swept set
  is wrong is worth more than any value drawn from it.
- **Open, with its price named:** $C_\lambda$ on $0<\lambda\le\lambda^*$ requires the degree-5
  fiber inequality over the **6-parameter** threshold family — three parameters wider than gate B,
  which is itself partial at budget.
- **Fifth float-versus-certified trap of the day, retracted inline:** a float Nelder–Mead run
  appeared to beat the paper's decimals by $+2.14\times10^{-11}$; certified re-evaluation at the
  optimiser's own digits is $2.11\times10^{-9}$ **below** them. Withdrawn per rule 5.

### `mm3/` gate C — a theorem, not just a negative: MONOMIAL TRANSFER

Frozen at `campaigns/2026-08-30T113838Z_c9532c07…_7e0e84591c56/` (pre-statement before compute,
`certified_landscape.json`, `invariance_check.json`, `audits.json`, `landscape.json.gz`, 23/23
checksums verified). **Outcome: no orientation in the swept set admits total $\le54$; certified
minimum $=55$. Nothing escalated, no record claim.**

- **MONOMIAL-TRANSFER THEOREM (the campaign's main finding).** For any signed permutation $g$ of the
  9 input coordinates, $C(F)=C(gF)$: applying $g$ to every wire maps inputs to $\pm$inputs (free),
  keeps $g(x\pm y)=g(x)\pm g(y)$ a single gate, and maps the 23 targets onto the 23 targets;
  $g^{-1}$ gives the converse. **Owner-verified structurally:** $d(\text{side})$ is the cardinality
  of a set of *sign-classes* of the 23 targets minus the 9 input directions, and $g$ maps
  $\{\pm v\}\mapsto\{\pm gv\}$ bijectively while fixing the input-direction set, so $d$ is
  invariant; full $C$-invariance follows because the whole synthesis problem is isomorphic under
  $g$ — which is what the identical DFS state counts measure ($33/116/66$ paper55, $132/81/672$
  mws59, $1152/1344/460$ stapleton60, $33/9/17$ sun56). **Consequence: the entire $48^3$ sandwich
  group is cost-invariant and the sweep collapses to 3 classes per decomposition.** This also
  retrospectively explains why session-2's order-3 orbit result was not hiding a larger landscape.
- **The sweep proved its own instruction redundant — and that is the result.** 1,658,880
  orientations (5 decompositions $\times$ $48^3=110{,}592$ signed-permutation triples $\times$ 3
  $\sigma$ rotations, 1160 s, checkpointed) produced **15 distinct data points**. The owner's
  instruction to sweep $48^3$ was therefore provably unnecessary — but the sweep is what
  *discovered* the invariance, which is now proven structurally rather than assumed. Both halves are
  recorded.
- **Owner-verified group counts.** Brute-forcing all $3^9=19{,}683$ sign patterns: ternary
  unimodular $3\times3$ matrices with $\lvert\det\rvert=1$ number **exactly 6960 = 48 monomial +
  6912 non-monomial**, confirming the agent's count. **Owner correction:** the agent wrote the
  unswept space as $6960^3\times3\approx1.0\times10^{11}$; the correct figure is
  $6960^3=3.372\times10^{11}$, hence $\approx\mathbf{1.0\times10^{12}}$ — an order of magnitude
  low, and it matters because that number is what justifies not sweeping.
- **Certified landscape (lower bounds, not achieved totals).** paper55 $=$ perminov58 (identical
  blocks, one object): all three $\sigma$ classes certified LB **55**. sun56: all three LB **55**.
  mws59: **56**. stapleton60: **58**. Note these are LBs, so the earlier session-2 *achieved*
  $55/57/57$ is consistent — paper55's $\sigma^1,\sigma^2$ classes now sit in $[55,57]$.
- **Three transcription defects caught BY the anchor requirement, all recorded inline per rule 5.**
  (i) A naive per-product $XUY^{-1}$ conjugation failed Brent — 12 failures on mws59 — replaced by
  the derived form $U'[r][(i_2,k_2)]=\sum (X^{-1})[i][i_2]\,Y[k_2][k]\,U[r][(i,k)]$; (ii) $\sigma^2$
  first coded as $(W^\top,U^\top,V)$ instead of $(W^\top,U,V^\top)$, 48 failures; (iii) Stapleton
  prints double-negated products ($M_0=(-t_3)\cdot(-B_7)$) and a draft encoded both signs
  separately, 2 of 729 cells failing. **Every one would have silently corrupted a plausible-looking
  landscape.** All five anchors (55, 56, 58, 59, 60) reproduced before the sweep, 729/729 each.
- **Counter is order-invariant by construction, not by a matching step** — $d$ depends on a *set* of
  sign-classes, so any summand permutation or per-summand sign flip leaves it fixed; additionally a
  pre-count machine check runs 8 pseudorandom relabelings per decomposition and aborts on any drift.
  All five passed, so relabelling cannot enter the landscape as variation.
- **Three new certified facts beyond the negative.** (a) **Sun's output stage is optimal**: sun56's
  $W$-factor floor is achievable at $d=16$ — the only achievable floor among all 45 side-instances
  swept — with an explicit verified 16-gate witness, so $C(W_{\rm fac,Sun})=16$ exactly and his
  $30=16+14$ is certified. (b) **stapleton60 certified LB 58** against a published 60, the first
  certified bound on that scheme here, from Appendix A data the agent **fetched and transcribed this
  session** (none existed on disk before). (c) mws59 certified LB 56 reproduced from an independent
  code path.
- **Rule 7, exactly as claimed:** no orientation in the swept set admits a certified total $\le54$
  and the minimum is 55. **Not claimed:** that 55 cannot be beaten. Unswept — the non-monomial
  ternary reorientations ($\approx10^{12}$), where monomial transfer does **not** apply and is
  therefore exactly where the landscape can move; other rank-23 decompositions; basis-change
  variants; $\mathbb F_2$-only schemes.
- **Strategic consequence for the gate's remainder.** The theorem discharges the monomial part *by
  proof* and simultaneously kills exhaustive enumeration for the rest — $10^{12}$ is a wall, not a
  sweep. The honest next instruments are structural (a transfer theorem for a non-monomial subgroup,
  or a bound on how much a non-monomial sandwich can help) or a **targeted search over sandwich
  entries** (the 27 entries of $(X,Y,Z)$ as ILP/SAT variables under unimodularity against the
  total-addition objective) — not a larger enumeration. Neither starts without a fresh
  pre-statement.
- **`2607.29291` closed in one line:** it yields no additional explicit $3\times3$ rank-23
  decomposition — its $\mathrm{GL}(3,2)^3$ + cyclic-trace + slot-matching machinery builds SAT
  witnesses for the Challenge CNFs — so nothing entered the sweep set, and the ring is recorded
  explicitly for all five members (all $\mathbb Z$).

### `mm3/` gate B encoding audit — the concern resolves, and the encoding is strengthened

The `2607.29291` warning (ten Heule–Kauers–Seidl formulas expected UNSAT are in fact satisfiable
because their CNFs *require* selected incidences without *forbidding* extras) prompted a bounded
written audit of `mm3/`'s floor-impossibility encoding. Frozen at
`campaigns/2026-08-30T113838Z_c9532c07…_7e0e84591c56/gate_b_encoding_audit.md`. **No
under-constraint; nothing escalated; $C(U)=13$, $C(V)=14$, $C(W\text{-factor})=14$ stand.**

- **The defect shape is INEXPRESSIBLE here, not merely absent** — a stronger answer than the
  inventory that was requested. `build_floor_cnf` has exactly **one** require family (coverage,
  $\bigvee_g Y[g,c]$) against **four** forbid families (at-most-one-class-per-slot;
  production-without-representation; representation-flag-without-production, i.e. full reification
  in both directions; representation-without-earlier-operand-availability). There are **no positive
  unit clauses at all**, only negative units forbidding unreachable classes — and the coverage
  require needs no companion forbid because slots are capped at one class with exactly $T=d(F)$
  slots for $d(F)$ classes, so "more than required" cannot be written down.
- **Both directions established, including the one Challenge-2 left unstated.** Model $\Rightarrow$
  circuit: read the schedule off the slots, operands available by the availability clauses. Circuit
  $\Rightarrow$ model: at $T=d(F)$ the floor lemma **forces** the gate-to-class map to be a
  bijection, so "every gate value is $\pm$ a needed target" is a **consequence, not a modelling
  assumption** — exactly where such an audit usually smuggles in the premise it is meant to check —
  and `prep()` enumerates all representation pairs so the matching $P$ variable is always settable.
  Hence UNSAT $\iff$ no $d(F)$-gate circuit.
- **Corroborated by a previously-measured behaviour, not merely consistent with one.** The encoding
  is deliberately too strong at $d+1$, which independently *predicts* the session-4 control already
  on record: raw $d{+}1$ UNSAT while aux-1 $d{+}1$ SAT and admitting the paper witnesses.
- **Tolerance question closed three ways, the third decisive.** The gate-C pipeline uses **no solver
  and no float** (exact Python ints plus `fmpz`, factor entries in $\{-1,0,1\}$, Brent sums bounded
  by 23). Gate B's HiGHS model has bounds exactly $[0,1]$, coefficients $\pm1$, RHS in $\{0,1\}$ —
  smallest nonzero magnitude $1.0$ against the $10^{-7}$ default, **seven orders of separation**, so
  the trap cannot arise. And the direction of any residual bug is safe by construction:
  `solve_floor_ilp` skips its `rhs=None` rows, which **relaxes** the model, and **an infeasible
  relaxation implies infeasibility of the tighter model**. Plus kissat DRAT replayed by `drat-trim`
  and `lrat-check` on the same instances.

### `delcap/` gate C — TARGET ACHIEVED: the published $d=1/2$ sandwich is improved at both ends

Frozen at `campaigns/2026-08-30T13:07:27Z_7988b619…_23445e8dbeeb/` — 40 files, 39/39 checksums
verified, `pre_statement.md` byte-identical to the pre-run version, and a freeze script that
**regenerates** its own anchor / second-route / support-guard evidence rather than transcribing it,
so those artifacts cannot drift from the claims. README appended $179\to458$ lines with the
pre-existing 179 verified byte-identical by `diff` after every edit.

- **The gate's stated target is hit.** A certified finite-blocklength improvement on
  Tavakoli–Nguyen–Bose's $d=1/2$ sandwich, **both ends, all 13 rows** ($q\in\{2,3,4\}$, $n\le8$).
  Truncation order **NONE** (full exact $q^n\times\sum_kq^k$ channel enumerated — no window, no
  alphabet cut), tail bound **exactly 0**, Arb 400 bits. Verdict every row
  `CERT_LOWER_BEATS_LBplus` **and** `CERT_UPPER_BEATS_UB`, by Arb-endpoint inequalities, margins
  $0.04$–$0.08$ and $0.08$–$0.29$ bits/symbol against widths $10^{-10}$–$7\times10^{-6}$: four to
  six orders of separation. `MACHINE-VERIFIED`.
- **Owner-verified at 300 bits:** $h_2(1/2)=1$ exactly; their closed forms
  $\mathrm{UB}=(1-d)\log_2q$ and $\mathrm{LB}_1=\mathrm{UB}-h_2$ reproduce their printed columns
  ($0.5$, $0.792481250361$, $1.0$ and $-0.5$, $-0.207518749639$, $0$ for $q=2,3,4$), confirming the
  agent's reading of Thm 1; and three spot-checked rows match the reported gains to the printed
  digit with the certified interval strictly inside $[\mathrm{LB}^+,\mathrm{UB}]$.
- **15 further rows at their own $(q,n,d)$: every certified interval strictly inside their sandwich,
  zero exclusions.** Plus **60 extension rows with no printed analogue** — 30 LO-CVB at
  $\delta\in\{1/20,1/2,4/5\}$ (ball radius $\le2.84\times10^{-120}$) and 30 Pinto–Ribeiro
  $C_{n,k}$ at $n=6..9$ — recorded as **first certified values, explicitly not comparisons**, since
  Morozov–Duman show those $\delta$ only as Fig. 2 curves.
- **The agent flagged its own undischarged requirement rather than letting it pass.** The
  $\mathrm{LB}^+/\mathrm{UB}$ columns are its **own exact-Arb recomputations** of TNB's Cor. 1 /
  Thm. 1 closed forms, **not** transcriptions of their printed Table I digits, so no row is compared
  at the digit level. Mitigated but not replaced: the recomputation is **not single-implementation**
  — route 1 is their closed form ($h_2$, $H_{\rm Bin}$, $\Delta_n$ from exact integer pattern
  counts), route 2 is the certified exact mutual information at **uniform input** from the full
  channel matrix, sharing none of that machinery, agreeing at $\lvert\text{diff}\rvert=0$ across 7
  points. **That agreement independently confirms TNB's own Cor. 1 claim that $\mathrm{LB}^+$ *is*
  the uniform-input rate** — a paper's structural claim verified as a by-product.
- **CONVENTION HYPOTHESIS FALSIFIED — by the very test the owner handed over. Owner hypotheses
  now 0-for-8.** All 18 TNB Table I entries transcribed first-hand. **TNB do NOT round up.** The
  single-offset fit is **infeasible under every convention** — round-half-up, ceiling and truncation
  bands all empty — and the reason is exactly the criterion derived for the MD case:

  | paper | printed $-$ computed | spread | ulp | single-$\Delta$ fit? |
  |---|---|---|---|---|
  | MD Table III (5 dp) | $[7.1001\!\times\!10^{-7},\,1.0399\!\times\!10^{-5}]$ | $9.6890\times10^{-6}$ | $10^{-5}$ | **exists** |
  | TNB Table I $\mathrm{LB}^+$ (3 dp) | $[-5.39\!\times\!10^{-4},\,+4.85\!\times\!10^{-4}]$ | $1.0240\times10^{-3}$ | $10^{-3}$ | **none** |

  One test, two papers, opposite verdicts, both owner-re-verified. TNB's signs are mixed nearly
  50/50 with every $\lvert\text{diff}\rvert\le5.39\times10^{-4}$, i.e. half an ulp — correct
  rounding applied **inconsistently row to row**, owner-confirmed on two rows that pin it: computed
  $0.824515\to$ printed $0.825$ is reachable *only* by round-to-nearest (truncation gives $0.824$),
  while computed $1.276539\to$ printed $1.276$ is reachable *only* by truncation (nearest gives
  $1.277$).
- **Consequence, taken honestly: the MD offset does NOT become typical and needs its own explanation
  after all.** Two of three papers round up (MD by measurement, PR by their own stated *"BA rate +
  tolerance"*), the third rounds to nearest. The agent declined to offer a replacement mechanism and
  explicitly noted its earlier $\sim10^{-6}$-relative note remains `INFERENCE` and is **not**
  strengthened by this test. The MD offset is therefore recorded as an **unexplained residual** with
  its direction and magnitude stated — not as a resolved convention.
- **What replaces the hypothesis is a better meta-finding, and it is the agent's.** Across the three
  papers this target consumes, printed tables differ from certified recomputation for **four
  distinct reasons**, none of which touches a theorem: MD Table III, conservative round-up (single
  $\Delta\le7.1\times10^{-7}$ plus 5-dp ceiling); PR, author-stated additive round-up; TNB Table I
  formula columns, inconsistent correct-rounding with no single $\Delta$; TNB's $C_{q,n}$ column, an
  unconverged BA estimate. **In finite-length deletion-channel work, printed tables are not reliable
  to their printed precision for several independent reasons — so certified recomputation is not
  redundant with reading the paper.** Worth more than the single-convention story that was hoped for.
- **Second escalation, re-verified and benign: TNB's $C_{q,n}$ column understates the finite-block
  capacity in 9 of 15 rows** by more than one 3-dp ulp (from $+1.307\times10^{-3}$ at
  $q{=}3,n{=}3,d{=}1/10$ to $+1.1097\times10^{-2}$ at $q{=}2,n{=}10,d{=}1/5$). Uses **only** the
  primal certificate — the exact mutual information of an explicit rational input in outward-rounded
  Arb — so it is independent of the dual. **No theorem is contradicted:** their stated ordering
  $\mathrm{LB}_1\le\mathrm{LB}_2\le\mathrm{LB}^+\le C_{q,n}\le\mathrm{UB}$ holds, and the
  certified intervals sit strictly inside their own $[\mathrm{LB}^+,\mathrm{UB}]$ in all 15 rows,
  which **confirms** the sandwich. **Owner-verified structure:** the excess is monotone increasing in
  $n$ at fixed $d$ ($q{=}2,d{=}1/5$: $+4.474\times10^{-3}$, $+7.298\times10^{-3}$,
  $+1.1097\times10^{-2}$ at $n=3,5,10$) **and** monotone increasing in $d$ at fixed $n$ in all three
  pairs tested. Both directions monotone is the unconverged-Blahut–Arimoto signature, and the
  decisive point is that **BA's primal rate climbs toward capacity from below**, so an unconverged
  run *understates* — our certified primal being a valid lower bound, exceeding their estimate is
  what an unconverged BA predicts. Recorded as: their $C$ column is an unconverged **estimate**, not
  a claim to printed precision (`INFERENCE`, from the monotone structure); **the positive statement
  is that this repository now holds tighter and higher certified lower bounds on finite-block
  capacity than the published BA column in 9 of 15 rows.** Framed as our bounds being better, never
  as their numbers being wrong.
- **Transcription authorized, with a sharper purpose than the pre-statement had.** Two of the three
  papers in this target are now known to print **conservative round-ups** rather than
  correctly-rounded values: Morozov–Duman's Table III (resolved above), and Pinto–Ribeiro who state
  it outright — *"reported value = BA rate + tolerance, an additive round-UP, so their printed
  numbers are valid upper bounds if their stopping criterion holds."* For an upper bound that is
  correct practice. TNB is the third and unknown, so transcribing their 18 printed entries now tests
  whether **round-up printing is a convention in this literature** rather than a quirk of one paper.
  Either outcome is informative: positive makes the MD resolution typical rather than anomalous;
  negative means the convention is not universal and the MD offset needs its own explanation.
- **Discipline applied without being told.** The channel-builder layer was anchored with *different*
  evidence from what it tests (TNB Examples 1–3 and the Thm 5 self-test: exact pattern histogram
  $\{1{:}4,2{:}2\}$, $\Phi_{1,2}(q)=1/q$, $\Delta_2(1/2)=1/4$) — rule 14, unprompted. Both gate-B
  anchors were re-certified through the **gate-C certificate code path** rather than the older
  pipeline, so the anchor exercises the layer that produced every gate-C number: $f(3,2)$ frozen
  $[1.4697819938,1.4697820261]\to$ gate-C $[1.4697819937559,1.4697819937572]$, contained and four
  orders tighter, additionally reproducing the independent frozen `rc28.py` dual at
  $\lvert\text{diff}\rvert=0$; $f(5,3)$ likewise contained.
- **A material deviation logged rather than quietly substituted, with a correct soundness
  argument.** The pre-statement fixed an mpmath-150dps locator; the rows used float64 because mpmath
  cost $\approx290$ s/row at $q{=}3,n{=}5$ and collapsed onto a near-degenerate $p$. Validity is
  genuinely untouched — the sandwich $[\text{primal},\text{dual}]$ contains the capacity for **any**
  candidate $p^*/D'$, so locator quality controls only tightness — and the achieved widths are
  *tighter* than the mpmath path's. Dropped points and the $m=10,15\to m=22$ substitution (their
  Table III actually prints $m=22$) likewise in `DEVIATIONS.md`.
- **One self-caught error worth the design note:** README anchor digits were first extrapolated past
  the 10 dp the manual run had printed; the freeze dry-run caught it and the README now carries the
  artifact's real 13-digit values. Precisely why the freeze script regenerates evidence instead of
  transcribing it.

### `delcap/` gate C — solver audit clean; an escalation that dissolved into a printing convention

- **Solver-tolerance audit: CLEAN, and re-verified rather than asserted.** `grep` for
  `scipy`/`linprog`/`highs`/`optimize` across all four gate-C scripts returns **0 hits** — no
  LP/QP/MILP anywhere in the certified path. Float appears in exactly two non-load-bearing places
  (a float64 numpy BA locating step; float64 ranking of $\Lambda$ subsets), both feeding an exact
  rational snap, after which every certified number is one outward-rounded Arb evaluation of logs
  of exact rationals at 400 bits. Certified per-symbol widths run $2.9\times10^{-10}$ ($n=2$) to
  $6.7\times10^{-6}$ ($n=8$) against sandwich margins of $0.02$–$0.13$ bits/symbol — four to six
  orders above any interval width, so no tolerance is ever compared against a quantity near its own
  size. One `SLSQP` call was an independent cross-check of the $n=3$, $d=1/2$ MI maximum only, never
  a source of a certified number; it agreed with the certified interval to $\sim9$ digits.
- **FIFTH instance of the session's failure shape — and the first one caught *before* it bit.** The
  agent audited the branch where the dual's reference distribution has $D'(y)=0$ while some
  $W(y\mid x)>0$: then $\mathrm{KL}=+\infty$, and silently **skipping** that term would yield a
  too-**small** "upper bound", i.e. an invalid certificate that looks fine. `cert_dual` returns
  $+\infty$ (no claim) instead, and full support is guaranteed two ways ($+1$ bump on any
  snapped-zero entry, plus an always-full-support uniform candidate, taking the min). **This was
  live, not hypothetical:** at $q=2$, $n=10$, $d=1/20$ the empty-output mass is $\sim10^{-13}$
  against a snap resolution $2^{-40}=9.094947\times10^{-13}$ — owner-verified, the resolution is
  $9.09\times$ **coarser** than the quantity it must represent, so the entry *would* have snapped to
  zero. Same shape as the HiGHS trap: a resolution coarser than the thing being resolved.
- **ESCALATION, frozen and unwritten: 19 of 36 Morozov–Duman Table III rows disagree.** Every
  disagreement is one-directional — our certified value **below** their printed value by
  $5.0\times10^{-6}$ to $1.04\times10^{-5}$, **above in zero rows**. Worst: $m{=}22,n{=}32$ printed
  $0.66391$ vs certified $[0.66389960,0.66389960]$; $m{=}5,n{=}1024$ printed $0.80279$ vs certified
  $[0.80277963,0.80277963]$. Interval widths $0.00$ at 400 bits (the quantity is $\log_2$ of an
  exact rational). The other **17** rows reproduce the printed 5 dp exactly, including all three
  $n\to\infty$ rows, which validates the Table I $E(m,w)$ transcription in aggregate. The agent
  verified its own evaluation by **two structurally independent routes** before flagging — Arb ball
  evaluation of summed rationals, and forming $L$ as one exact big rational then $\log_2$ in mpmath
  at 80 dps — identical to all printed digits.
- **OWNER RULING: the conclusion stands, the explanation is FALSIFIED.** The agent's `INFERENCE`
  was "last-digit rounding in a printed table". Tested in exact decimal arithmetic: $0.66389960$
  rounds to $0.66390$ under **both** half-up and ceiling, and $0.80277963$ to $0.80278$ — neither
  reproduces the printed $0.66391$/$0.80279$, and both differences ($1.040\times10^{-5}$,
  $1.037\times10^{-5}$) **exceed one ulp** at 5 dp. A printed-rounding artifact cannot exceed one
  ulp, so the explanation is withdrawn per rule 5.
- **What the arithmetic actually indicates.** Those two differences differ by only
  $3.0\times10^{-8}$ across rows as dissimilar as $(m{=}22,n{=}32)$ and $(m{=}5,n{=}1024)$ — the
  signature of a **systematic additive term**, not rounding and not noise. But 17 rows match
  exactly, including every $n\to\infty$ row, so it is not a global constant either: it is a
  **finite-$n$ term, present in some rows and absent in others.**
- **Not a refutation, and no claim is contradicted.** A larger printed value is the *weaker* upper
  bound, so their converse stands. The $19$–$0$ direction split is itself decisive against random
  error (two-sided $p=3.8\times10^{-6}$ under a fair coin). Framed correctly this is a **result**:
  in those 19 rows our certified values are strictly tighter than the printed table, with
  zero-width intervals — the first certified enclosure of those rows, improving published
  precision, not a complaint about the paper. **The 19 rows stay out of the index until the
  separator is found.**
- **RESOLVED, and it is benign — the escalation dissolves into a printing convention.** The
  classification was ordered before any mechanism, and it dismantled the agent's split *and* the
  owner's hypothesis.
  * **The 19/17 split carried zero information.** Sorted by gap, the rows partition *perfectly* at
    the sort boundary: the split was a cut at $5\times10^{-6}$ — a half-ulp window the agent chose
    itself — through a **continuum**. The two classes interleave on every property named ($m=5,22,23$
    on both sides, $n$ from 1 to $\infty$ on both, $\lvert\Lambda\rvert$ from 5 to 24 on both);
    gap correlations are noise at $n=36$ ($n$: $+0.40$, $mn$: $+0.28$, full-set: $+0.20$, $m$:
    $-0.08$, $\lvert\Lambda\rvert$: $-0.03$).
  * **Direction count corrected against the agent's own earlier report: 36–0, not 19–0.** Every row
    has printed $>$ certified, minimum gap $+7.1\times10^{-7}$.
  * **The owner's "systematic additive term" read is FALSIFIED — hypotheses now 0-for-7.** It was
    inferred from the two gaps that happened to be quoted sitting $3.0\times10^{-8}$ apart; those
    are merely the top two of a continuum (third $1.0238\times10^{-5}$; **35 distinct values among
    36 rows** at $10^{-8}$ resolution). A top-of-distribution proximity artifact — the owner
    generalized from the two extreme data points in a report instead of asking for the distribution.
  * **The whole thing reduces to one inequality: the spread is less than one ulp.** Gaps span
    $[7.1001\times10^{-7},\,1.0399\times10^{-5}]$, spread $9.68899\times10^{-6}<10^{-5}$. Under
    ceiling printing, gap $\in[\Delta,\Delta+\mathrm{ulp})$, so the **minimum** gap caps $\Delta$
    above and the **maximum** gap bounds it below: $\Delta\in(3.990\times10^{-7},\,
    7.1001\times10^{-7}]$, owner-re-derived and matching the agent's band to the digit. **Exactly
    determined by the extremes, not fitted** — and non-empty *precisely because* spread $<$ ulp.
  * **Evidence against a formula difference, not merely absence of evidence for one:** a formula
    term would be $n$- or $m$-dependent and would **break** the single-$\Delta$ fit; the fit holds
    across all 36 rows. The residual $\le7.1\times10^{-7}$ absolute at rate values $0.55$–$0.82$ is
    $\approx0.9$–$1.3\times10^{-6}$ relative — the size of a 6-to-7-significant-digit intermediate.
  * **SSOT drift found in this repo's own index, same class as the `physics/` one found this
    morning.** `RESULTS.md` carried **6 target rows for 7 targets** — `delcap/` had **no row at
    all**, while the summary paragraph above the table asserted that "`delcap/` has a certified
    Blahut–Arimoto pipeline plus two documented negative findings". Prose claiming a state the
    index lacks is exactly the `physics/na-compiler` pattern from earlier today, recurring here
    within hours of being diagnosed there. Row added with its full evidence, and the drift itself
    recorded in the row so the omission is not silently repaired. **Standing lesson: an index that
    is asserted complete in prose must be counted, not read.**
  * **Final framing: NOT a discrepancy.** Morozov–Duman's values agree with ours to within
    $7.1\times10^{-7}$ absolute in all 36 rows, and they printed with **round-up** — the correct,
    conservative convention for an upper bound, since a rounded-up upper bound stays valid. The
    apparent $10^{-5}$ disagreement is 5-dp ceiling rounding atop a $\sim10^{-6}$ relative precision
    difference. Nothing is contradicted, their converse stands, and the convention reflects well on
    them. Recorded as **36 rows independently recomputed and certified** (zero-width at 400 bits,
    `MACHINE-VERIFIED`; the convention identification `INFERENCE`) — the first certified enclosures
    of those rows. Both retractions kept inline per rule 5.

---

## 2026-08-30 — `kg/` closed: the $3\times$ excess was our defect, the paper stands

The longest-running open thread in the repo is resolved, against the owner's own
repeatedly-tested hypotheses, and the resolution is a repo-side bug rather than a
published error.

- **Root cause: plain derivatives where the paper requires the Euler operator.** Paper
  line 857 defines $D:=t\,\mathrm d/\mathrm dt$; line 1816 gives the recursion
  $\Phi_{k+1}:=D_t\Phi_k-t\,\partial_x\partial_y\Phi_k$; line 1826 states
  $\rho_j=D_t^j\rho$. Since $D_t^jt^m=m^jt^m$, the closed form is
  $\rho_j(t)=(t-3^js_3^2t^3+5^js_5^2t^5)/V$. The code used $\partial_t^j\rho$, evaluated
  at $t=+1$ where $t=-1$ was required — both of the agent's reported values at
  $\theta=\pi$ matched plain derivatives at $t=+1$ to eleven digits. Error sizes:
  $\rho_2$ was $27.92\times$ too large, $\rho_3$ was $3.376\times$ too small and
  sign-flipped. The recursion was never underdetermined; it was being read at lines
  1826–1854 where the c-table is *printed* rather than at 1816 where it is *defined*.
- **Corrected sweep lands inside the paper's bound.** avg$(S)=123.377609$ against
  $\le126.80385221$, ratio $0.97297998$ — $2.7\%$ slack, exactly what a certified upper
  bound looks like. $B_3=14.2459830 < 14.4424366$. **Hypothesis (iii) "the paper's
  certified constant errs" is DEAD, not merely unadvanced.**
- **Headline strengthens.** Corrected slack to break-even rises from the retracted
  $2.1576124\times$ to $3.7851739\times$, and our own route now lands below the paper's
  boundary rather than merely consuming its cited tail — so the `CITED-DEPENDENCY` tail
  is independently corroborated, at `COMPUTATIONAL-EVIDENCE` grade only.
  $K_G\le1.7818413238804372880$ unchanged.
- **Independent structural confirmation, from the paper's own averaging step.** Every
  exponent of $\rho_j$ is odd, so all $\rho_j$ are odd, so $S$ is $\pi$-periodic; line
  1861's replacement of $(1/2\pi)\int_0^{2\pi}$ by $(1/\pi)\int_0^\pi$ is valid *only*
  under that periodicity. Plain derivatives have mixed parity and cannot satisfy it.
  Measured: Euler gives $S(0)=S(\pi)=171.6916$; the defected run gave $141.649$ vs
  $1711.417$, a $12\times$ symmetry violation previously logged as an unexplained curiosity.
- **A subagent refutation of the paper was caught and voided.** A claimed certified lower
  bound $\lVert D^3H\rVert^2\ge290.382 > 208.584$ was re-verified by the owner and found
  to sum Parseval past the grid's exactness edge: the grid is the odd-parity half of the
  *square* $a,b\le251$, so $b_m$ is exact only for $m\le251$ and the missing rows carry
  mixed sign. Reproducing the procedure one truncation level down turned $0.777362899$
  into $8755.04$. Retained from it: $\lVert D^3H\rVert^2\ge0.77736289857822393$.
- **Owner calibration.** Structural hypotheses ended 0-for-4 (numerology twice, then an
  $L^2(\mu)$ factor-3 route falsified by a clean 14/14 unit test), plus two
  instrument-specification errors: conflating the paper's intermediate with its
  conclusion, and calling the Parseval route "the decider" when line 1861's inequality
  gives it zero power over avg$(S)$. What actually worked was process — the
  second-route requirement, the escalate-before-concluding rule, the audit-non-transfer
  ruling, and owner re-verification of every agent number, which caught the void
  refutation and an incorrect $B_3=14.2202$ before either reached the index.
- **All three closures landed, and every figure was owner-recomputed from the frozen
  `leg2euler_*.npy` arrays rather than taken from the agent's report**: the Euler c-table is
  committed inside `d3h_sweep49_asm.py` and a clean-process rerun reproduces
  avg$(S)=123.377609$ ($25\to49$ increment $10^{-6}$, spectral); a fresh independent
  mpmath-40-dps evaluator puts assembly-vs-jet at $6.6\times10^{-15}$, retracting the
  agent's $1.2\times10^{-3}$ probe-truncation figure; and $\pi$-periodicity is a standing
  sweep assertion at $3.96\times10^{-15}$. Owner also derived and verified a **stronger**
  check the agent never ran: $S$ is even as well as $\pi$-periodic, which forces the full
  reflection $S(\theta)=S(\pi-\theta)$ — 24 independent equalities across the grid, not
  one — measured at $7.936\times10^{-15}$, with $\min S=1.266545$ falling exactly at the
  symmetry point $\theta=\pi/2$.
- **Label held.** The row stays `COMPUTATIONAL-EVIDENCE`: the $\theta$-quadrature is
  uncertified float, so spectral self-convergence and $10^{-15}$ structural agreement do
  **not** earn `MACHINE-VERIFIED`. Promotion needs certified widths through the quadrature.

---

## 2026-08-30 — all seven targets carry committed evidence; three firsts, three dissolved false alarms

Closing state of the launch session. Every target now has a frozen campaign
artifact and a verdict; no target is left at protocol stage.

- **`mceliece/` gate A: PASS, Apon's mechanism confirmed by measurement.** Seven
  synthetic binary-Goppa instances $(m,n,t)$ — $(6,64,3)$, $(6,64,4)$,
  $(6,64,5)$, $(7,128,3)$, $(7,128,6)$, $(8,256,3)$, $(8,256,5)$. On every
  instance and every held count $c=1\ldots(2t+4)$: $N_{\rm fam}(c)=2t+3-c$
  **exactly**, reaching $0$ precisely at $c=2t+3$; every per-(point, order
  $j<4$) family-restricted block has rank **exactly 1**. Zero refutation events,
  zero degenerate held points, $\Delta_{p,q}\ne0$ throughout. 7 seconds on one
  core. The force of this comes from construction: the flag rows are built
  explicitly from Apon's Definition-4 map $T(A)=Q_A+\rho_A Z F$ (eqs. 17–20), so
  Lemma 7's "at most one condition per held point on the $(2t+3)$-family" is
  **derived by measurement, not assumed**. Scope stated: a finite theorem over
  those instances; it does **not** establish $c_{\rm need}>2t+3$ at
  `mceliece8192128` or any cryptographic scale. Frozen at
  `mceliece/campaigns/2026-08-30T01-42Z_9D0FC9E7/`, with the agent's two
  self-found algebra bugs (Rabin Frobenius exponent; GRS multiplier $\lambda$ vs
  $\lambda'$) and the VOIDED conclusion they had produced preserved in
  `PROBE_RESULTS.md` rather than deleted.
- **`omega/` gate B rung 2 certified, and a reusable second-order finding.** At
  VXXZ24's released point (6759 params): $\omega\le2.3715538358544617$, a valid
  rigorous bound that sits $2.0258545\times10^{-6}$ **above** their published
  $2.37155181$ and so does not reproduce their six-decimal rounding. The
  shortfall is entirely the ln-space dual residual of the **shipped** Lagrange
  multipliers (worst $3.0767\times10^{-6}$) against interval widths of
  $1.08\times10^{-19}$ — the arithmetic is five orders of magnitude tighter than
  the gap. Decisive: with tight multipliers the published value **is**
  recoverable, margin $3.7942377\times10^{-9}$. So this is a statement about
  released certificate data, **not** about the validity of the bound and **not**
  about any error in the paper. **The positive result leads:** the two-rung
  ladder order survives rigorous arithmetic — rung 1 certifies to
  $2.3713400836689$, rung 2 to $2.3715538358544617$, gap
  $2.1375219\times10^{-4}$ — the first two-rung rigorous comparison in the
  combination-loss line, and both endpoints also sit below DWZ23's $2.371866$.
  **General finding, first measurement of its kind:** the *price of rigour* on
  this method's artifacts is $1.084\times10^{-6}$ (Alman25) and
  $2.026\times10^{-6}$ (VXXZ24), in both cases dominated by shipped-multiplier
  looseness rather than arithmetic. A cheap multiplier re-centring probe was
  attempted and honestly reported as insufficient (residual scattered across
  shapes, not a uniform offset), with the full fix named as a follow-up campaign.
- **`delcap/` delivered the first rigorous enclosure table for the
  Fertonani–Duman genie-aided quantities.** All 13 published Table II 2-dp
  entries enclosed in exact interval arithmetic, and every one sits **above**
  its certified interval — so they are valid upper bounds but not tight, with
  the round-up-on-an-unconverged-dual mechanism reproducing all 13. One **new**
  certified row not in the published table: $f(11,5)\in[2.3888120200,\,
  2.3888120204]$. Rubinstein–Con's $0.3745$ at $d=0.68$ needs $n=28$ and is
  honestly reported unreachable here ($U_{\rm cert}(28,0.68)=0.13168$ against
  the row's $0.1199$; band certified to $n=16$), rather than presented as
  confirmed at reduced $n$. Durable asset: a certified-Arb Blahut–Arimoto with
  rational snap and Arb dual certificate, duality gap $\le10^{-11}$ for
  $n\le7$ — machinery this literature did not have.
- **`oct-rank/` closed all three gates.** Gate B's rank-24 hunt is a
  well-powered clean negative over 273 probes, and the border-rank
  instrumentation is what makes it trustworthy: a continuation curve showed
  residual holding at $\approx1.3\times10^{-15}$ while norms grew
  $13.7\to1135.7$, but on **25-column** rows, i.e. inter-term re-weighting
  rather than border-rank drift; norms snapped back at $s=0$ and the true
  24-column endpoint ran away. An $s=0.9$ run reaching $1.747\times10^{-15}$ was
  **self-escalated and rejected** by the agent as re-convergence to the rank-25
  certificate. Gate C PASS: $\mathrm R_{\mathbb R}(\tau)=7$ re-verified both
  directions, and $\frac52n-2$ confirmed sharp at $n=2$ (rank 3) and $n=4$
  (rank 8). No claim that $\mathrm R_{\mathbb R}(T_{\mathbb O})>24$.
- **`mm3/` closed all three gates and its ladder.** Gate B now rests on **three
  independent decision procedures** (DFS, HiGHS infeasible, CaDiCaL UNSAT) plus
  witnesses. Gate C: no $\le54$ in the valid $\sigma$-orbit; landscape
  $\{55,57,57\}$, the first machine-checked per-orientation optimality map for
  rank-23 $3\times3$; pure $(U,V)$ swaps excluded **by proof** ($B\cdot A\ne
  A\cdot B$), not by budget. Full predecessor ladder re-verified over
  $\mathbb Z$ in one convention: Stapleton-60, MWS-59, Perminov-58, Sun-56, and
  the 55-paper.
- **Three apparent problems in the literature were raised and all three
  dissolved.** Sun-56 versus the 55-paper was a *different rank-23
  decomposition*, not a contradiction (zero common product triples). The
  59-versus-57 flag was the **agent's own tally slip**, found by machine recount
  and retracted inline — convention checked first and matched, record history
  $60\to59\to58\to56\to55$ unchanged, no underclaim by any author. The FD10
  Table II gap was a round-up convention on an unconverged dual. **No refutation
  of any published claim is asserted anywhere in this repository.**
- **Interval-arithmetic hygiene became repo knowledge, from two independent
  discoveries.** `Omega` root-caused that python-flint's `arb(lo, hi)`
  constructs (midpoint, **radius**) rather than an interval — silent, and it had
  produced a phantom zero-width collapse. `Delcap` then found the complementary
  trap: a genuine radius of $5.8\times10^{-89}$ on a value near 1 reads as
  exactly $0.0$ if width is computed by subtracting float64-converted endpoints,
  so the owner's own "assert nonzero width" diagnostic was itself degenerate and
  had to be corrected to a **precision sweep**. Both were broadcast; all
  ball-arithmetic targets audited clean afterwards, including a precision sweep
  on the repo headline (widths $6.10\times10^{-121}$ → $7.08\times10^{-30}$ as
  precision drops $512\to96$, margin positive throughout).
- **Owner corrections to the record this wave.** The `kg/` row previously said
  an independent computation put $b_1$ at $0.881598$, making the paper's value a
  conservative lower bound; that was **wrong** — the figure was a partial
  evaluation missing the $A_{0,1}$ term, and the correct closed form matches the
  paper's printed $b_1$ to $4.86\times10^{-17}$. The headline stands; its
  mechanism is now correctly attributed entirely to the authors' conservative
  choice of $\gamma$. Also corrected: `omega/`'s double-count (endpoint improved
  $1.0836\times10^{-6}$, slack term $16{,}424.7\times$ tighter) and its own
  "100× / margin unchanged" summary (actual: margin improved to
  $2.1173\times10^{-4}$).
- **Successor agents dispatched** for the two highest-value open items:
  `OmegaRung2` (completed, above) and `KgD3H` — an independent re-derivation of
  $\lVert D^3H\rVert$, which would lift the headline's last CITED-DEPENDENCY and
  clear its $2.7\times$ width qualification in one step. Two earlier attempts at
  that task yielded without producing anything on disk and were replaced.

## 2026-08-30 — headline: a certified improvement on the Grothendieck constant; `mm3/` complete through gate C

- **HEADLINE RESULT.** `kg/` gate A certifies
  $K_G\le1.7818413238804372880\ldots$, which is **stronger than the bound
  arXiv:2608.11158 states** ($1.7818666069360661$) by $2.528\times10^{-5}$ —
  an improvement over Krivine of $3.72654\times10^{-4}$ where the paper claims
  $\ge3.4737\times10^{-4}$. Mechanism: the paper's sufficiency condition is
  $\gamma+\Delta_H<b_1$, so the largest admissible value is
  $\gamma^\*=b_1-\Delta_H$; the authors reported a round, conservative $\gamma$
  instead of the maximal one. Nothing in their mathematics is wrong; the last
  $2.5\times10^{-5}$ was simply never extracted.
  Certified at 256 bits outward-rounded:
  $b_1^{\rm low}=0.88157382204959954858188766422240828$ (half-width
  $2.46\times10^{-36}$), $\Delta_H^{\rm upper}=1.59045454374832\times10^{-5}$
  (head $1.13288599276897\times10^{-5}$ + tail $4.57569\times10^{-6}$),
  $\mathrm{margin}_{\rm low}=1.2508504162\times10^{-5}>0$,
  $\gamma^\*=0.881557917504162$. Frozen at
  `kg/campaigns/20260830T010702Z_2a373482_385dbaab2fd3/`.
  **Owner re-check (independent, Arb 300 bits, this session):** confirmed. The
  decisive check is that $\pi/(2\gamma_{\rm paper})=1.7818666069360660912$
  reproduces the authors' own `certificate.json` value $1.7818666069360661$
  **exactly**, which pins the identification of their $\gamma$ and shows
  $\gamma^\*$ is real unclaimed slack rather than an arithmetic artifact.
  **Caveats recorded as load-bearing, not footnotes:** the verdict is "PASS
  with tail input CITED" — $\lVert D^3H\rVert\le14.44243664663976457$ and the
  A-grid `grid_M251.json` are imported CITED-DEPENDENCY from
  `github.com/trishullab/grothendieck-bounds`, so this is an arithmetic-level
  improvement on the authors' own scheme and coefficients, not an independent
  certification of the scheme; and the dominant enclosure width is that
  imported tail at $4.58\times10^{-6}$, only $2.7\times$ below the margin
  rather than the $10\times$ the owner required, kept as `width_qualification`
  in the artifact. Gate C was reprioritised accordingly: re-derive
  $\lVert D^3H\rVert$ **first**, which would lift the last cited dependency and
  clear the width rule in one step, and only then attempt the septic extension.
  Escalate-not-adjust locked in: if the re-derivation disagrees with 14.4425 the
  agent reports the discrepancy and stops.
  Also recorded: $K_{\rm DR}=1.6769\ldots$ is `[REPORTED]` after a full
  dead-link audit (`bound2.dvi` → 301 to cse.umn.edu, zero status-200 CDX
  captures; Davie 1984 unobtainable) — the retrieval failure the gate ladder
  predicted, mirroring what stalled `../math/ns/` on NRS and ESS. No gate
  depends on it, verified.
- **`mm3/` complete through gate C, and an apparent inter-paper contradiction
  dissolved.** Gate B's optimality verdict is now backed by **three independent
  decision procedures** — exhaustive floor DFS, HiGHS ILP infeasible, and
  CaDiCaL 1.5.3 UNSAT (396/1440, 468/1807, 481/1820 vars/clauses at
  $T=d(F)$) — plus the paper's own circuits as $d(F)+1$ witnesses, giving
  $C(U)=13$, $C(V)=14$, $C(W\text{-factor})=14$, output $=28$. The agent
  **retracted its own over-reading** in which aux-1-impossibility was thought
  part of the certificate; the correct obligation is floor-impossibility plus
  witness, and aux-1 at $d+1$ is possible as the witnesses require.
  Gate C: **no $\le54$ total exists in the valid orientation orbit.** Swept set
  stated exactly — the $\sigma$-orbit $\{\mathrm{Id},\sigma,\sigma^2\}$,
  $\sigma:(U,V,W)\mapsto(V,W^\top,U^\top)$ of order 3; pure $(U,V)$ swaps were
  excluded **by proof, not budget**, since they compute $B\cdot A$ rather than
  $A\cdot B$. Certified totals $\mathrm{Id}\to55$, $\sigma\to57$,
  $\sigma^2\to57$, so the paper's orientation is the certified optimum of its
  own orbit and the landscape $\{55,57,57\}$ is the first machine-checked
  per-orientation optimality map for rank-23 $3\times3$. Complete over the
  valid orbit, partial over abstract orientation space — so "no 54" is not a
  universal claim (rule 7).
  The Sun tension resolved as a **decomposition mismatch, not a contradiction**:
  his 13-addition circuit is real (verified twice, from his shipped `verify.py`
  and his printed §4, 729 identities, 0 failures) but computes his own $U$-map,
  and his decomposition is not a reorientation of `cr58_cn122`-as-printed —
  factor-column multiset equality fails, and still fails up to per-column signs,
  under the $T$ involution, and entry-transposed, with **zero** common product
  triples. Bonus: the same procedure independently confirms Sun's own lower
  bounds, including his reported `V_aux1_at_12_possible: False` over all 338
  auxiliary directions. Two published papers corroborated by one procedure.
- **`rs-pe3d/`: the first exact 3-dimensional product-expansion ratios ever
  computed.** Nine instances, each $\delta$ exact over $\mathbb F_q$; per-witness
  exactness then *proven* for all nine — every line-support tuple of cost below
  $\delta$ verified absent and at least one at $\delta$ verified present
  (`optimality.json`). Values: $q=13$, $s=(2,2,4)$ → $1/2,\,3/8,\,1/4$ at
  $t=(1,1,1),(1,1,2),(1,1,4)$; $q=31$ → $(2,3,5)\,3/5$, $(2,3,10)\,3/5$,
  $(2,6,5)\,1$, $(3,10,2)\,3/5$; $q=61$ → $(2,3,5)\,3/5$, $(3,4,5)\,1$.
  Labelled MACHINE-VERIFIED **finite-window**: best-of-window with the window
  (point supports of weight $\le3$) stated, global $\rho$ OPEN for heavier
  supports. First balance data: at $q=61$, $\beta=1.67\to\rho=1$ versus
  $\beta=2.5\to\rho=3/5$ — visible degradation with imbalance, explicitly not a
  collapse claim. That is the axis where the proved 2-dimensional theorem
  carries a balance hypothesis Conjecture 4.2 drops.
- **`omega/` gate B stage (a) reproduced; stage (b) failed honestly and was
  reopened.** The Alman25 rung's released verifier (~1.3k LOC MATLAB, 11 files)
  transcribed to Python with a 24855/24855 group-count match on the first
  attempt; $\omega$ in file $2.3713389005434182$ against published $2.371339$
  ($9.95\times10^{-8}$, display rounding); max inequality violation
  $1.137\times10^{-10}$, equality $2.390\times10^{-11}$, an order inside their
  own refine tolerance; the Schönhage constraint violated by $-6.2\times10^{-15}$
  of half-ulp float noise. Stage (b)'s interval aggregation produced an
  arithmetically impossible $R_{\rm sum}=-4.27$ where the true value is
  $\approx+2.8$, and the agent correctly labelled that **transcription-plumbing
  evidence, not a statement about the bound**, claiming no endpoint. Reopened
  with a per-block containment test (each Arb enclosure must contain the stage-(a)
  float) to localise before any margin is computed. Live lead: the shipped
  Lagrange multipliers leave a Lemma-1 residual at order $3.7\times10^{-6}$, so
  the eventual honest statement is $\omega\le2.371339+\varepsilon$ with
  $\varepsilon$ dominated by that residual, not by float noise — and such an
  endpoint would still sit well below VXXZ24's $2.37155181$.
- **Process items.** A new SSOT rule was added to `README.md` after agent `Mm3`
  replaced `mm3/README.md` with a campaign summary, deleting the gate ladder:
  target READMEs are **append-only below the contract**, agents add a
  `## Current state` section and never delete the contract. The ladder was
  restored with the agent's content retained verbatim beneath it. Corrections
  now carry attribution in both directions — four owner-side fixes
  (`omega/`'s "float-only" line, `rs-pe3d/`'s $12\mid q-1$ condition,
  `mceliece/`'s infeasible ladder, `kg/`'s runtime estimate, the last replaced
  by the agent's measured 545 s $b_m$ loop and ~10 min end-to-end figure) and
  one agent-side fix (`rs-pe3d/` had misattributed its own false "Anchor 1"
  premise to the owner's README; the README contains no such claim).
- **`mceliece/` still has no frozen verdict**, and was redirected: run the
  waterfall census at **small $t$ first**, where $2t+3$ is single-digit and the
  sweep is seconds, rather than waiting on the $m=10..12$ ladder. Its probe
  work already earned three results — a Rabin exponent bug ($Z^{2^d}$ instead
  of $Z^{q^d}$) that would have silently voided the squarefreeness guard; the
  GRS multiplier correction to $\lambda'=G(a_i)^2/\Pi'(a_i)$, with the earlier
  probe conclusion explicitly VOIDED; and the analytic finding that $F$ is
  unique (evaluation injective since $D=n-2t-1<n$) so $F'\equiv0$ is impossible
  for a nonzero row, which makes Apon's $\Delta_{p,q}\ne0$ a **genuine**
  genericity condition to be measured per instance rather than a formality.

## 2026-08-30 — first gate verdicts; two repo-level findings

All seven owning agents ran. Every one committed `pre_statement.md` before its
first computation, as required. Verdicts and findings, newest work first:

- **`oct-rank/` gate A: PASS.** arXiv:2608.16649 independently re-verified end
  to end, campaign frozen at
  `oct-rank/campaigns/2026-08-30T00:04:39ZZ_gateA_A6A4B0CB/`. Exact-rational
  re-check of the rank-25 Krawczyk certificate: $K=0.31955674<1$,
  $\mathrm{off}_{\max}=1.314\mathrm e{-14}$, worst per-equation margin
  $+6.804\mathrm e{-7}>0$, **all 512 per-equation margins strictly positive**;
  outward-rounded Arb re-evaluation strict. Independent radius re-derivation
  [DERIVED]: hypotheses hold for all $\rho<\rho^\*\approx3.129\mathrm e{-6}$
  (60-iteration bisection; contraction-limited, all three flags flipping
  together just above), so the paper's $\rho=10^{-6}$ carries $\approx3.1\times$
  slack. Lean 4 replay: `leanprover/lean4:v4.29.1`, mathlib `5e932f9`,
  `lake build` 2398 jobs clean, axiom audit of all key theorems exactly
  `[propext, Classical.choice, Quot.sound]`, no `sorry`. Verdict
  MACHINE-VERIFIED `[REPRODUCED]`; the mathematics remains
  CITED-DEPENDENCY. Gate B (rank-24 attack) running.
- **`mm3/` gate A: PASS.** arXiv:2607.28676 reproduced exactly; campaign frozen
  at `mm3/campaigns/2026-08-29T235739Z_c3a32a44_7c3ed00afa49/`. All 729 Brent
  identities over $\mathbb Z$ in exact `fmpz`: 27 unit, 702 zero, 0 failures;
  counts $13/14/28=55$ confirmed; the paper's SLP, printed factor blocks and
  expanded products agree exactly (three independent presentations). Perminov's
  `cr58_cn122` (commit `98ba522`) decoded with his own loader gives a tensor
  identical to the paper's after the stated C-permutation, 0 mismatches.
  MACHINE-VERIFIED `[REPRODUCED]`. Gate B (ILP + SAT/PB optimality decision)
  running.
- **`omega/` gate A finding: the current $\omega$ world record is not
  independently checkable from the published paper.** As of 2026-08-29
  arXiv:2608.16884 gives the reformulated program (Eq. 11, $q=5$,
  $\ell^\*=4$, ~7M parameters, Lemma 1 $H^{\max}$ scheme) but publishes **no
  parameter table and no ancillary files**, stating only that the authors *"are
  preparing a repository in which we will release the verification code and our
  discovered solution."* Audited negative, re-runnable: no repository from any
  of the ten authors, `google-deepmind/alphaevolve_results` holds only the
  May-2025 notebook, text search for `2.371177` finds nothing. **This is a
  statement about artifact availability on a date, not a claim the bound is
  wrong.** Framing correction adopted: an earlier line in `omega/README.md`
  saying the method's published constants are "all float-only" was **wrong** and
  has been replaced — all three predecessor rungs released parameters on OSF
  (Duan–Wu–Zhou `osf.io/dta6p` $2.371866$; Vassilevska Williams–Xu–Xu–Zhou
  `osf.io/7wgh2` $2.37155181$, 6759-dim float64 vector; Alman et al.
  `osf.io/mw5ak` $2.371339$, 24855-dim float64 vector), so they are
  float-reproducible but never interval-enclosed. Gate B redirected there, under
  a mandatory two-stage transcription control (reproduce the published float
  value first, only then switch the same code path to Arb) because the released
  verifiers are MATLAB and a transcription bug could masquerade as a failed
  enclosure.
- **Repo-level finding — the checkability spectrum.** Same field, same month,
  both AI-assisted, opposite ends: arXiv:2608.11158 ($K_G$) ships
  `github.com/trishullab/grothendieck-bounds` with frozen Arb certificates,
  per-coefficient mid/rad pairs for $b_3..b_{251}$, an explicit
  $\lVert D^3H\rVert$ certificate and a stated margin $1.2554\mathrm e{-5}$;
  arXiv:2608.16884 ($\omega$) ships a promise. This is direct evidence for the
  repo's positioning thesis recorded yesterday — AI-assisted results are
  arriving faster than they are being certified, and the variance in artifact
  quality is enormous.
- **`delcap/` — a false positive correctly dissolved, and a reduction correctly
  falsified.** Two process results worth more than either row:
  1. An apparent discrepancy in Fertonani–Duman's Table II (their $f(3,2)=1.48$
     vs a certified capacity $1.46978\ldots$) was **not** an error. Measured
     across 12 tabulated pairs the offset is *irregular*, $+0.0006$ to
     $+0.0127$, never a uniform $+0.01$, and all 12 published 2-dp values sit
     above the certified intervals. Mechanism identified and reproducing all
     12: their stated round-up rule applied to a **not-fully-converged
     Blahut–Arimoto dual**, which approaches capacity from above. So the Table
     II entries are valid upper bounds, not capacity estimates. **No literature
     correction warranted**; recorded in the campaign artifact as a calibration
     datapoint, deliberately not in any README.
  2. A proposed symmetry reduction ($S_n$-equivariance $\Rightarrow$ capacity of
     $W_{n,k}$ equals that of the type channel on $n+1$ inputs) was **tested
     before being built on, and falsified**: type-channel capacity is strictly
     below full-channel capacity at every tested $(n,k)$ except a coincidental
     agreement at $(3,1)$ — e.g. $(3,2)$: $1.4698$ vs $1.2539$; $(4,3)$:
     $2.1699$ vs $1.5850$. Reason: the $S_n$ action on inputs induces a
     non-permutation action on outputs, because deletion preserves order and the
     output is a subsequence. Reduction discarded; no $n=28$ claim rests on it.
     Path to Rubinstein–Con's $0.3745$ is now their own published space
     reduction, which is their stated contribution.
  Also: Rubinstein–Con's $0.1221$ lower bound reproduced from their released
  `BDC_Lower_Bounds` pipeline ($0.12213624912193806$ here vs $0.12208469671934222$
  in their notebook vs $0.1221$ in the paper), labelled `[REPRODUCED]` /
  COMPUTATIONAL-EVIDENCE only, because their `k_probs` clipping at $10^{-300}$
  makes the entropy term non-rigorous. Certified-Arb BA pipeline now closes to
  gap $\le10^{-11}$ for $n\le7$; $n=8$ costs ~100 s at gap $10^{-4}$; $n=28$
  unreachable densely.
- **`kg/` gate A in flight, with the error budget now the binding constraint.**
  Parameters are print-precise in the paper (§4 Eq. 7:
  $\eta=0.136419125$, $s_3=0.34101124$, $s_5=0.05276111$,
  $\vartheta=\eta/\sqrt V$, $V=1+s_3^2+s_5^2$; $\gamma_{\rm paper}=0.881545409$),
  handled as exact rationals. Because the margin is only
  $1.2554\mathrm e{-5}$, the verdict was tightened to
  $\mathrm{margin}_{\rm low}=b_1^{\rm low}-(\gamma+\Delta_H^{\rm high})>0$ with
  all enclosure widths reported beside it, and the composite Gauss–Legendre
  quadrature error is being bounded **explicitly** per panel rather than left
  implicit in ball arithmetic. Scope declared in advance: the tail input
  $\lVert D^3H\rVert\le14.4425$ is imported CITED-DEPENDENCY from the authors'
  frozen certificate, the $N=251$ tail is re-derived here, and the verdict is
  labelled "PASS with tail input CITED" on every restatement; full
  re-derivation is scheduled in gate C where the machinery is needed anyway.
- **`mceliece/` gate A imminent.** All four dispute papers plus the upstream
  distinguisher read in full. The flag definition is pinned to Apon's Eq. (3)
  and Step-3 system (6), which Vedenev's Eq. (13) matches; Saarinen's three
  reductions are relation-generation-layer and are correctly excluded from the
  Step-3 assembly, with his Eqs. (29)/(30) reading of Apon used as the
  cross-check anchor. Convention risk closed by adopting Apon's own §2
  experiment protocol, so neither disputant can dismiss the result as a
  mismatch. The measurement is **two** family-restricted quantities —
  $N_{\rm fam}(c)=\dim\{A\in E[Z]_{\le2t+2}:\text{flag conditions hold }\forall
  a\in S_c\}$, predicted $2t+3-c$, and the per-(point, order) rank of the flag
  restriction on the family, predicted exactly $1$ — and ambient rank is never
  measured, so the spurious-refutation mode is designed out. Instance guards
  machine-checked, including Apon's differential identity
  $\Pi F'+\Pi'F=\kappa G^2F^{(2)}$ coefficient-wise.
- **`rs-pe3d/` gate A starting**, pre-statement committed with Definition 4.1,
  Conjecture 4.2, Theorem 2.1/4.11 and the AI claim quoted verbatim, and the
  exact quantity fixed as $\rho_{\rm inst}=\min \mathrm{wt}(M)/\delta(M)$.
  Instance-selection correction adopted from the agent: this repo's README had
  said "choose $q$ with $12\mid q-1$", which is an imprecise proxy; the real
  requirement is that $q-1$ admit **three pairwise coprime divisors $>1$**
  ($q=31$ qualifies via $30=2\cdot3\cdot5$ despite $12\nmid30$, and
  $\mathbb F_{25}$ is a valid field but $|\mathbb F_{25}^\star|=24$ admits no
  such triple, so $q=25$ is correctly dropped). Adopted census
  $q\in\{31,61,109,181,241,421,601\}$.
- **Owner steering issued this session** (all recorded in the agents' inboxes):
  blocked `delcap/`'s unproved symmetry reduction before it could carry a
  result; required `omega/` to state its availability finding with a date and
  the authors' stated intent quoted, and imposed the two-stage transcription
  control; warned `oct-rank/` that a CP fit with diverging factor norms
  certifies **border** rank, not rank, and required joint residual/factor-norm
  reporting plus a Krawczyk failure-cause split; required `kg/` to report
  enclosure widths beside a $1.2554\mathrm e{-5}$ margin and to bound quadrature
  error explicitly; required `mceliece/` to declare its convention and to
  measure family-restricted rank only.

## 2026-08-29 — seven targets opened; seven attack agents launched

- Scan complete. All six research slices reported; distilled, owner-checked
  output is [`docs/CS_FRONTIER_SCAN_2026-08-29.md`](docs/CS_FRONTIER_SCAN_2026-08-29.md),
  raw scout output preserved verbatim and non-authoritatively in
  `docs/scan-raw/`. Aggregate listing coverage recorded there: arXiv `cs.CC`
  2026-04…08 full, `cs.DS` 2026-07 full, `cs.IT` 2026-07 (414) + 2026-08 (362),
  ECCC TR26-001…154 full index, IACR ePrint 2026 (~id 1740–1823).
- **Convergence:** two independent scouts plus the owner's own `cs.CC` sweep
  arrived at the same top candidate, the 2026 Grothendieck-constant band.
  `LearnComm` ranked it #1 (`krd-cert/`), `MetaComplexity` referred it out
  rather than develop it, and the owner had already found it in the
  2026-08 listing. Opened as `kg/`.
- **Owner first-hand reads this session** (V-owner, not delegated):
  arXiv:2608.11158 (abstract, §1, §2 in full), arXiv:2608.16649,
  arXiv:2608.16884, arXiv:2607.28676, arXiv:2305.07156, ePrint 2026/1747,
  2026/1810, 2026/1786, the ECCC latest-reports feed (TR26-145…154), and
  arXiv `cs.CC` 2026-08 listing entries 1–36 and 51–82.
- **First machine-verified fact.** Arb at `ctx.prec=400`: $6\pi/11$ and
  $\pi/(2\log(1+\sqrt2))-10^{-4}$ both lie strictly inside $(1.7,1.8)$, so the
  tenths digit of $K_G$ is $7$; interval width
  $0.06851798532420916352206252\pm3.37\mathrm e{-27}$; and $6\pi/11>1.67696$,
  beating the cited Davie–Reeds constant. This certifies the **corollary** of
  arXiv:2608.11158, not its bounds. Recorded in `RESULTS.md` as the headline
  with that limitation stated.
- **Field shift recorded, owner-verified from the ECCC pages:** AI-discovered
  proofs are now credited inside primary TCS venues (TR26-146 credits GPT-5.6
  Sol and Claude Fable 5 for the upslice scheme and its proof; TR26-150 credits
  ChatGPT 5.6 Sol with crossing the Polishchuk–Spielman sum-of-rates barrier in
  2D; arXiv:2608.11158 and arXiv:2608.16884 both attribute their results partly
  to engineered AI systems). Strategic consequence adopted for this repo: AI
  results are arriving faster than anyone certifies them, so the repo positions
  as the **certification layer**, not as a competitor in generation. TR26-150's
  own higher-dimensional conjecture ships with an AI proof the authors state is
  unverified; the $\omega$ record ships with no rigorous enclosure. Both became
  targets.
- **Seven targets opened**, each with a README carrying a numbered gate ladder,
  pre-campaign primary-read requirements, disjointness statement, and numeric
  pass/fail thresholds: `mceliece/`, `kg/`, `omega/`, `delcap/`, `oct-rank/`,
  `rs-pe3d/`, `mm3/`. Ranked shortlist and rationale in the scan document.
- **Toolchain ready before launch:** `cs/.venv` (Python 3.14.3) with
  python-flint 0.9.0 (Arb + exact rationals + finite fields), sympy 1.14.0,
  numpy 2.5.2, scipy 1.18.1, python-sat 1.9.dev15, highspy, z3-solver 5.1.0,
  mpmath 1.3.0. Smoke test passed: Arb 1/3 at 200 bits, exact `fmpq(6,11)`,
  pysat UNSAT detection, HiGHS import.
- **Seven owning agents launched in parallel** (`Mceliece`, `Kg`, `Omega`,
  `Delcap`, `OctRank`, `RsPe3d`, `Mm3`), each under a strict single-folder
  ownership contract: an agent may edit only its own target folder; this
  session (Main) remains the single writer of `README.md`, `RESULTS.md`,
  `PROGRESS.md` and `docs/`. Cross-folder code sharing is forbidden — a
  duplicated small helper beats a dependency between concurrent writers.
  Mandated order of work: primary reads → `pre_statement.md` committed →
  implement → run → freeze `campaigns/` artifact. Refutations of published
  claims must escalate to Main before being written anywhere as established.
- **Corrections issued during the scan** (detail in the scan document §
  "Corrections issued"): (i) `MetaComplexity`'s proposed `oct-rank/`
  lower-bound gate by integer-witness enumeration is **unsound** — real tensor
  rank is semialgebraic, so no enumeration bounds it below; the target's gates
  were rewritten to make the upper side witness-producing and the lower side an
  audit, with "a failed search yields no claim" stated explicitly. (ii)
  `ProofSAT` failed (exit 1) mid-correction after detecting its own material
  misattribution — the ~200 TB DRAT / 68 GB / 800-core facts belong to Boolean
  Pythagorean Triples (Heule–Kullmann–Marek 2016), **not** van der Waerden
  $W(2,6)$, and Kouril–Paul (2008) shipped no certificate; a dedicated
  `ProofSATFix` pass was dispatched, and only `rs-pe3d/` was promoted from that
  slice pending it. (iii) `LearnComm` misstated the Krivine sufficiency
  condition as $\gamma+\sum_{m\ge3}|b_m|\gamma^m<b_1$; the correct condition
  (owner-read) is $\gamma+\Delta_H<b_1$ with $\Delta_H=\sum_{m\ge3}|b_m|$, the
  $\gamma$-powered sum being the majorant $M(\gamma)$, a different object.
  `kg/README.md` carries the correct form and the launch brief warns the agent
  not to reintroduce the error. (iv) The scout's `mceliece/` instance ladder
  ($\lambda\in\{9,10,11\}$, $t\in\{72,80,88\}$) is **infeasible** — at $m=9$,
  $t=72$ we get $mt=648>512\ge n$, so $k=n-mt<0$ and no such Goppa code exists;
  replaced with a nine-row ladder that satisfies $k>0$ throughout while
  preserving the intended $2t+3$ threshold range.
- Rejected with reasons (scan document § "Rejected this round"): meta-complexity
  structural questions, pseudorandomness constructions, transformer/LLM
  expressivity, the Williams $\sqrt{t\log t}$ and Cook–Mertz constant audits
  (attractive but no published number to refute), and quantum-adjacent
  complexity (reserved to `../physics/`).

## 2026-08-29 — repository opened; six-scout TCS frontier scan launched

- `cs/` created as the third sibling attack repository beside `../math/` and
  `../physics/`, scoped to **theoretical computer science and information
  theory**. Same four-layer SSOT contract; `RESULTS.md` authoritative index,
  this ledger chronological, `<target>/README.md` per-target state,
  `<target>/campaigns/` immutable artifacts.
- Discipline: `../math/README.md` rules 1–7 inherited verbatim, plus six
  TCS-specific rules (model statement; finite-$n$ is never asymptotic; measured
  not assumed complexity; searches carry certificates; conditional results keep
  their hypothesis; `[REPRODUCED]` vs `[DERIVED]`). Recorded in
  [`README.md`](README.md#discipline).
- Selection bias declared up front: prefer targets whose progress is a
  machine-checkable finite certificate on one workstation over targets whose
  only currency is asymptotics or floating-point benchmark. Rationale: a solo
  repo loses scaling races and can own the certification niche.
- Disjointness constraints fixed before the scan: quantum error correction,
  quantum hardware, and quantum compilation are reserved to `../physics/`
  (`qldpc-dec/`, `msd/`, `shadows/`, `na-compiler/`, `fss-bb/`, `qlops/`) and
  `../math/qec/`; the ten `../math/` targets are likewise off-limits.
- Six research agents launched in parallel over disjoint field slices:
  `MetaComplexity` (meta-complexity / space-vs-time / circuit + algebraic
  lower bounds), `FineGrained` (fine-grained complexity, $\omega$, graph
  algorithms), `InfoCoding` (capacities, coding theory, LDC/LCC, exact LP
  bounds), `ProofSAT` (proof complexity + certified combinatorial search),
  `CryptoHardness` (assumptions, concrete-security audits, proof systems),
  `LearnComm` (learning theory, communication/query complexity, LOCAL).
  Each was required to page through real listing indexes (arXiv cs.CC / cs.IT /
  cs.DS / cs.DM / cs.CR / cs.LO, ECCC recent reports, IACR ePrint 2026) and to
  report the listing URLs it opened, so scan coverage is auditable.
- Citation rule imposed on every agent: an identifier may be printed only after
  its landing page was fetched and its exact title quoted; otherwise the
  citation is dropped or tagged `UNVERIFIABLE`. Raw output lands in
  `docs/scan-raw/` and is **non-authoritative** — scout tags are
  subagent-verified, not owner-verified. Any claim entering a formal gate
  statement gets an owner first-hand primary read first.
- Writer discipline for this session: this session (Main) is the single writer
  of `README.md`, `RESULTS.md`, `PROGRESS.md`, and the distilled scan document;
  agents write only their own file.
