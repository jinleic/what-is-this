# N2 VERDICT.md — upward Gröbner/elimination campaign (OctRankNext)

Campaign: `campaigns/20260901T115801Z_8d6683fb_39c752e7eadd/` — gate
`n2_groebner_infeasibility`. Prereg committed BEFORE any build/compute:
commit `028f5880ddde433653ef035dc44828039604193e`, prereg sha256
`46e596bc80b91c290694bb639c1d50f4bda2a99050c31243b0a7353a3cf423a4`,
byte-identical copy in this dir as `pre_statement.md` (PROVENANCE.md).

## Headline

**BUDGET-FAIL — quantified, first-class resource negative. The
Gröbner/elimination route to rank-13 infeasibility is measured to be out of
reach: the SMALLEST pre-registered calibration system (K-2, the
KNOWN-INFEASIBLE tau r=6 system: 66 variables / 48 equations — a factor 4
smaller than the 273-variable main system) blew its 30-minute S1 cap while
still inside F4 degree 9, at an 8,080,476 x 74,499,764 Macaulay matrix and
22.65 GB resident memory. Per the pre-registered S1 rule ("if (b) fails to
certify infeasibility within cap: BUDGET-FAIL for the whole campaign — the
instrument cannot scale even to the known case; no claims"), the main S2
attempt was NOT run. rank((L_1, L_i, L_j)) stays OPEN in {13, 14}.**

**Absence of a 13-witness is NOT evidence for 14; absence of an
impossibility argument is NOT evidence for 13.**

The published window 18 <= R_R(T_O) <= 25 is untouched; no published work is
refuted. This campaign makes NO rank claim: a resource negative is a
statement about the ROUTE, not about the tensor.

## What was fixed before compute

Everything load-bearing: the two system formulations (E: gauge-extended
273 vars / 218 eqs; E1: E + pinned term, 273 vars / 220 eqs), the sound
gauge-elimination construction and its semantics, the three primes
(32003, 65537, 2147483647) and their order, the monomial order (DRL /
grevlex, msolve default), the char-0 protocol (rank-14 claimable only via
Q-mode {1} or the stated clause-(ii) strength), stage caps (S1 30 min CPU
each, S2 6 h/prime, 24 h campaign), the four controls, and the verdict
mapping. No parameter was adjusted after seeing results.

## Instrument (built inside the campaign, post-commit)

- msolve **0.8.0**, built from the source tarball with sha256
  `319ba0de67dca967dea40cb6e4dacf44eab387c2f0c2416b71842c11affbadfd`
  (byte-identical to the artifact Sage 10.7 pins), linked against homebrew
  GMP + MPFR and a locally built FLINT 3.6.0.
- Built binary sha256: `.libs/msolve` =
  `1e8a6725f0eef8d0c158f432eaad1e4bb4539205f50827048612fa778225bc21`;
  wrapper `msolve` =
  `2e2944a0eec3415d70ff3530b50e4fbe74240a0d97e9c93284ff0d813568742a`.
  (Build tree moved out of the run dir after the runs to keep the frozen
  artifact lean; BUILD_LOG + these hashes + configure.out/make.out remain.)
- Instrument smoke test: the Sage spkg-configure vector
  (`x-y, x*y-1` over Q) reproduces the expected msolve output byte-for-byte;
  a 2-variable GB smoke (`x^2+y^2-1, x*y` mod 101) returns the correct
  3-element reduced basis in 0.18 s. The tool works; the SYSTEMS are the
  obstruction.

## Systems serialized (exact integers; hashes frozen in n2_systems.sha256)

| file | role | vars | eqs |
|---|---|---|---|
| s1a_tau_r7.txt | K-1 accept-known-feasible | 77 | 48 |
| s1b_tau_r6.txt | K-2 reject-known-infeasible (calibration target) | 66 | 48 |
| s1c_tau_r6_corrupted.txt | K-3 one-coordinate perturbation smoke | 66 | 48 |
| k4a_tau_r7_ext.txt | K-4a gauge-extension feasibility | 91 | 62 |
| k4b_tau_r6_ext.txt | K-4b gauge-extension infeasibility | 78 | 60 |
| s2e_tf_r13_ext.txt | **S2 main E** | 273 | 218 |
| s2e1_tf_r13_ext_pin.txt | **S2 main E1** | 273 | 220 |

Prime-field variants (char line rewritten, same generators) exist for
32003 / 65537 / 2147483647 as declared. Generation-time anchors passed:
tau has 12 nonzeros (3 slices x 4), TF has exactly 24 = 2x12 (blockdiag
doubling), all entries in {-1, 0, 1}.

Size correction (disclosed, rule 5): the prereg wrote "E1: 275 variables";
the pin uses the EXISTING a0_0 / c0_0 variables, so E1 is 273 vars / 220
eqs. Recorded in code with a comment; no threshold or domain moved.

## Measured costs — the negative, quantified (in-process/OS measurements)

K-2 (tau r=6, mod 32003, `msolve -g 1`, i.e. leading monomials only — the
CHEAPEST possible question), exact command:

    ./msolve -v 2 -g 1 -f s1b_tau_r6_p32003.txt -o k2_lm_try3.txt

F4 round-by-round trace (real | cpu seconds), from `k2_try3_capped.log`:

| deg | pairs selected | matrix | new | zero | time (real \| cpu) |
|---|---|---|---|---|---|
| 4 | 192 | 384 x 2123 | 192 | 0 | 0.02 \| 0.02 |
| 5 | 2684 | 2688 x 18981 | 616 | 364 | 0.20 \| 0.20 |
| 6 | 16996 | 24120 x 159817 | 2112 | 2976 | 5.85 \| 4.88 |
| 7 | 93826 | 167308 x 1243634 | 7604 | 17312 | 145.82 \| 121.55 |
| 8 | 519761 | 1171307 x 9566028 | 30908 | 79536 | 3802.54 \| 3593.23 |
| 9 | 2900755 | **8080476 x 74499764** | (in progress at cap) | | killed |

- Stage reached: F4 **degree 9**, matrix 8.08e6 x 7.45e7 (~6e14 entries in
  the dense sense; 0.00% density), still in linear algebra when the cap
  bound.
- Memory high-water observed: **22.65 GB** RSS at the kill point in this
  attempt; **35.93 GB** RSS observed in the earlier attempt of the same
  system (both from OS `ps` snapshots, recorded in
  k2_final_ps_snapshot.txt — labelled OS-derived, used ONLY as the resource
  measurement this negative is about, never as a mathematical number).
- Degree-8 alone consumed 3802 s wall / 3593 s CPU — i.e. **more than twice
  the entire 30-minute S1 cap in a single F4 round**, on a system four times
  smaller than the main target.
- Three independent attempts of K-2 all failed to produce output within
  caps (0-byte outputs: k2_lm.txt, k2_lm_try2.txt, k2_lm_try3.txt);
  wall times 90 s (probe), 1860 s, 1740 s. K-1 (tau r=7, 77 vars) had
  already exceeded a 120 s probe window.
- S2 (E / E1 at 273 vars) was **NOT run**: the pre-registered S1 rule
  forbids it once K-2 fails, and the measured degree-9 blow-up at 66
  variables makes the 273-variable case unreachable by many orders of
  magnitude within any pre-registered budget.

## Controls status (honest, incomplete by cap)

| control | required | observed |
|---|---|---|
| K-1 accept known-feasible (tau r=7) | must NOT return {1} | INCOMPLETE — exceeded its probe window; no output |
| K-2 reject known-infeasible (tau r=6) | must certify {1} within 30 min | **FAILED BY CAP** (degree 9, 22.65 GB, no output) — this is the campaign's decisive measurement |
| K-3 one-coordinate perturbation smoke | run to completion, outcome recorded | NOT RUN (S1 rule fired first) |
| K-4 gauge-soundness (extended tau) | feasibility/infeasibility survive extension | NOT RUN (S1 rule fired first) |

No control ever returned a FALSE infeasibility, so no quarantine condition
was triggered; equally, no control validated the pipeline's positive
direction. Therefore NOTHING in this campaign is reported as mathematical
evidence about ranks — only as measured resource facts about the route.

## Reproducibility of this negative (assignment item 3)

Everything needed is frozen here: `n2_systems.py` (deterministic generator,
exact integer arithmetic, generation-time anchors), `n2_systems.sha256`
(hashes of all seven serialized systems + a meta sidecar per system),
`n2_systems_manifest.json` (var/eq counts), the exact msolve commands in
this verdict and in the logs, `BUILD_LOG` + `configure.out` + `make.out`
(build recipe), the binary and tarball sha256s, the full F4 traces
(`k2_verbose.log`, `k2_try3_capped.log`, `k2_try2.log`), the OS memory
snapshot (`k2_final_ps_snapshot.txt`), and the 0-byte output files that
witness non-termination. A successor can re-run the identical commands
against the identical inputs and will hit the same wall.

## Verdict mapping (deviation disclosed)

The prereg mapped BUDGET-FAIL to `FROZEN-INCONCLUSIVE`. The campaign is
closed as **FROZEN-NEGATIVE** on the owner's standing instruction
(2026-09-01, Main): a cap-bound outcome with quantified measurements is to
be frozen as a first-class negative ABOUT THE ROUTE. The semantic content is
exactly the pre-registered BUDGET-FAIL: no rank claim, no impossibility
claim, {13,14} unchanged. This wording change is recorded here rather than
silently applied.

## Rule-7 scope sentence

Covered: msolve 0.8.0 (pinned build) on the seven declared integer systems
at the declared primes/order, with the S1 calibration cap binding on the
tau r=6 control; measurements as tabulated. NOT covered, hence not decided
by anything here: the rank of (L_1, L_i, L_j) (still OPEN in {13,14});
whether the E/E1 systems are feasible or infeasible; any other solver,
order, prime, or reformulation (e.g. signature-based/F5, sparse/structured
elimination, SAT/SMT over R, Positivstellensatz certificates); the frozen
floor 13 and upper 14 (both intact, untouched); border rank; complex
decompositions; R_R(T_O) and its published window.

## Named next action (for the owner, not claimed here)

The elimination route needs a fundamentally cheaper formulation before it is
worth any further budget: exploiting the blockdiag(tau,tau) symmetry to cut
variables (the 26 gauge unknowns are dwarfed by the 247 CP unknowns), or a
structured/sparse solver, or attacking a proper sub-family whose
infeasibility would still imply rank >= 14. The measured degree-9 wall at 66
variables is the calibration number any such proposal must beat.
