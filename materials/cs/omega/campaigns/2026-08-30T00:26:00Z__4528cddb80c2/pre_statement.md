# Pre-statement — `omega/` campaign (fixed before any run)

Written 2026-08-30T00:03Z, after primary reads, before any computation in this campaign.
This file fixes the exact quantities, precision, rounding modes, budgets, certificate
format and pass/fail thresholds. Thresholds may not be adjusted after results are seen.
Any re-scoping is recorded in a numbered addendum below with its reason; the original
statement stays visible.

## Primary reads (completed before this statement)

All four papers read first-hand in full; excerpts quoted in
`campaigns/` artifacts. Cited by exact arXiv ID + version, titles verified on
landing pages this session:

- **[REC] arXiv:2608.16884v1** — Dupont, Eisenberger, Kozlovskii, Mehrabian, Ruiz,
  See, Zhou, Alman, Vassilevska Williams, Balog, *"Improving the matrix
  multiplication exponent with modern optimization and AlphaEvolve"* (2026-08-17).
  TeX e-print downloaded and grepped; no ancillary files; assets/ contains only logos.
- **[PRED1] arXiv:2404.16349 (v3, 2026-08-19)** — Alman, Duan, Williams, Xu, Xu,
  Zhou, *"More Asymmetry Yields Faster Matrix Multiplication"* (SODA'25).
  TeX e-print: §7 numerical results; released code+params at
  `https://osf.io/mw5ak/` (view-only link printed in the TeX).
- **[PRED2] arXiv:2307.07970 (v1)** — Vassilevska Williams, Xu, Xu, Zhou, *"New
  Bounds for Matrix Multiplication: from Alpha to Omega"* (SODA'24). TeX §8;
  released code+params at `https://osf.io/7wgh2/` (view-only link in TeX).
- **[PRED3] arXiv:2210.10173 (v3)** — Duan, Wu, Zhou, *"Faster Matrix Multiplication
  via Asymmetric Hashing"* (FOCS'23). TeX §7 has footnote: "Data and code are
  available at https://osf.io/dta6p/". Downloaded `fmm-code.zip`.

## Exact object of gate A (the record bound)

The record constructs a feasible point of the program [REC, Eq. (11)]:

```
minimize  Ω
subject to  E_total + M_total · Ω ≥ 2^{ℓ*−1} · log₂(q+2),
            all free variables lie in the domains of [REC, §2.2].
```

with hyperparameters q = 5, ℓ\* = 4 (both stated in [REC §3, first sentence]).
Any feasible point certifies ω ≤ Ω [REC, Theorem 1, quoting Alman et al. 2025].

**Discovery (this session, before Computing):** the optimizing parameters are not
published anywhere [REC] can point to. Verified facts:
- [REC] contains no parameter tables (TeX grepped).
- arXiv e-print has no ancillary files.
- [REC §4] verbatim: *"We are preparing a repository in which we will release the
  verification code and our discovered solution."* No such public repository exists
  as of 2026-08-30 (author GitHub accounts, google-deepmind/alphaevolve_results,
  GitHub/OSF general search all checked this session).
- All three predecessors DID release: [PRED1] `osf.io/mw5ak` →
  `code_matrix_mult.zip` with `data/W1.00_2.371339.mat` (flat float64 vector,
  24 855 params + MATLAB verify harness); [PRED2] `osf.io/7wgh2` → `rmmcode.zip`
  with `data/K100_2.37155181.mat` (6 759 float64 params + verify harness);
  [PRED3] `osf.io/dta6p` → `fmm-code.zip` with
  `data/power8_dup_2.371866.mat` (MATLAB MCOS objects; params embedded).

**Gate A finding is therefore the dated artifact-state statement** (per Main's
directive, 2026-08-29): *as of 2026-08-29, the parameters attaining
ω < 2.371177 are not published at any precision and no release has appeared, so
the bound is not independently checkable from the paper as published.* This is a
statement about the current state of the artifact, with a date. It is NOT a
claim that the bound is wrong, NOT a claim of impropriety, and NOT a claim that
it is unverifiable in principle. The authors' stated intent, verbatim from
[REC §4]: *"We are preparing a repository in which we will release the
verification code and our discovered solution."*

Search ledger (negative result, auditable and re-runnable when the release
lands; all searches 2026-08-29):
1. arXiv e-print Tarball of 2608.16884 — no ancillary files; `assets/` contains
   only 3 DeepMind/Google logo PNGs; TeX grepped for parameter tables: none.
2. GitHub accounts of reachable co-authors checked via the GitHub API:
   `orbitingflea` (Renfei Zhou; 14 repos, none omega-related; homepage feed has
   no post about a release), `EmilienDupont` (8 repos; none), `marvineisenberger`
   (not found), `bkozl` (not found), `emilien-dupont` (1 unrelated repo).
3. `google-deepmind/alphaevolve_results` — contains only the May-2025
   `mathematical_results.ipynb` for the AlphaEvolve technical report; no omega
   optimization artifacts.
4. GitHub full-text search for "2.371177" and repo search for
   `omega matrix multiplication exponent` and for description text: 0 hits.
5. OSF search for the record: no project.
6. Web search for the exact bound string with verification/repository keywords:
   no release page; secondary coverage (aiweekly, emergentmind, officechai)
   repeats the paper text only.
I will NOT re-run an optimizer and present its
  ℓ\*=3, q=5, 4th power of CW₅.
- [PRED2] ω ≤ 2.371552 point, file `K100_2.37155181.mat`, 6 759 float64 params,
  ℓ\*=3, q=5.

Each parameter vector v is an exact rational point after reading its float64 bits
exactly (float64 = dyadic rational). The claimed bound Ω(v) is the value of
transcribed [PRED1]/[PRED2] verification programs. The enclosure is computed over
the interval hull conv(v) — the unit box around the float vector obtained by
width-2⁻⁵² outward rounding of each entry (candidates: 2⁻⁵² and 2⁻⁴⁰; chose
2⁻⁵², narrowed from 2⁻⁴⁰ only if the enclosure is too coarse; that is a precision
budget, not a threshold change).

Precision and rounding mode (fixed):
- Base arithmetic: `python-flint` 0.9.0 `arb` ball arithmetic, mid-point precision
  200 bits, with every final comparison made on **outward-rounded** endpoints:
  upper bounds are rounded UP, lower bounds rounded DOWN, at 100-bit endpoint
  precision, and every inequality consumed is an interval inequality.
- Entropies and logs: all logs base 2; interval `log` from `arb`, outward as above.
- Max-entropies: every inner maximization bounded rigorously FROM ABOVE by a
  Lagrange/dual certificate of the [REC §2.5 Lemma 1] type: for a max-entropy
  target with marginals ρ_W, a dual vector y with y(a)>0, same marginals, and
  λ₀+λ_X(a_X)+λ_Y(a_Y)+λ_Z(a_Z) attaining |log y(a) − g(a)| ≤ ε gives
  H^max(ρ) ≤ H(y) + 2ε (Lemma 1, machine-verified in our reproduction too). When
  the released data file stores its own `max_entr_dist` / `dist_max` components,
  those serve as the y's and the stored Lagrange multipliers serve as the λ's;
  every residual is computed in interval arithmetic outward.
- Monotone/affine structure: where a consumed expression is affine in the float
  params (most of [PRED1]/[PRED2] programs: sums of H terms with fixed weights),
  bounding H on the box conv(v) by concavity of H gives a rigorous upper bound via
  tangent-plane/dual bounds at the float point — [REP verified in code].

**Search budgets (fixed):** the box conv(v) is NOT sampled; any subdivision used
for the dual-residual bounding has explicit budget ≤ 10⁶ leaf boxes and a
stop condition (all leaves resolved to ≤ target width 2⁻⁴⁰ or certified above
by concavity tangent bounds). No randomized search. No daemon loops. One
low-priority process at a time (`nice -n 10`), BLAS/OMP threads pinned to 1.

**Gate B pass/fail thresholds (fixed now):**
- PER RUNG (separately for 2.371339 and 2.371552): PASS = certified upper endpoint
  of the resulting Ω enclosure < published rounding, i.e. < 2.371339 (rung 1),
  < 2.371552 (rung 2), computed from the released parameters. This is the
  *reproduction* test: **does the released parameter file really certify the
  published number, when ALL arithmetic is interval and every inner max is bounded
  rigorously above?**
- MARGIN-reporting rule: if PASS but the certified Ω exceeds the published decimal
  while still < published rounding, the reported margin is the certified gap to the
  published rounding target. If the enclosure FAILS to close below the published
  number, that is a genuine finding: the released verify harness (which allows max
  violation ≤ 10⁻⁶ or 10⁻⁹) accepted a point whose rigorous enclosure is wider than
  the published 6-decimal claim; the exact step where the float-to-interval gap
  opens must be reported (which constraint, at what interval width).

**Certificate format (fixed):** each gate-B run produces a JSON file with
- run id, UTC timestamp, code hash (git-style sha256 of src),
- exact input vector digest (sha256 of float64 bytes),
- per-constraint list of interval widths and certified endpoints,
- the maximum residual constraint name,
- the final certified Ω interval endpoints (decimal strings, ≥ 12 significant
  digits),
- every dual certificate (λ's, y's, ε's) consumed,
or an explicit `FAIL` field naming the first constraint where the interval gap
opened and the certified upper endpoint that resulted.

**Never claimed:** that 2.371177 (or either predecessor rung) is *re-proved
from scratch* here. Published bounds are CITED-DEPENDENCY; our enclosure outputs
MACHINE-VERIFIED interval statements about *the released parameter vectors*,
never about the unreleased record.

## Addendum A1 (2026-08-30T00:25Z, per Main directive — framing correction)

The omega/README.md gate-B phrase "the method's published constants are all
float-only" is INACCURATE and Main is correcting it on their side. Accurate
state: **all three predecessors released float64 parameter vectors plus
verification code on OSF; they are reproducible in floating point, but were
never enclosed in interval arithmetic.** Every artifact and report produced by
this campaign uses the accurate version. Gate B's purpose is therefore unchanged
and sharper: it asks whether the released float64 points *rigorously* certify
the published 6-decimal numbers once every inner maximization is bounded above
rigorously and every evaluation is outward-rounded — the first rigorous
enclosure of any rung of the combination-loss ladder, if it passes.

## Addendum A2 (2026-08-30T00:25Z, per Main directive — mandatory two-stage
transcription protocol for gate B)

Transcription of ~2.4k LOC of MATLAB to Python under interval semantics is the
primary threat to gate B. A failed enclosure that is a transcription bug, reported
as "the record does not certify", is the worst output this repo could emit.
Mandatory protocol, before any interval result is trusted:

(a) **Float64 reproduction stage.** Transcribe the target rung's MATLAB
verification program to plain float64 Python, with no interval logic, and
reproduce the published bound value from the released parameter vector to within
the paper's own stated tolerance. The released harness ships a 1e-9 refine mode
and warns at 1.1e-6, so the target tolerance is |reproduced − published| ≤ 1e-9
(|r − p| ≤ 1e-6 counts as weak match and is reported as such). Record the
reproduced value next to the published one. Only a stage-(a) match is evidence
that my transcription is faithful.

(b) **Interval stage.** Only then switch the SAME code path to Arb with outward
rounding (change of arithmetic backend, not of program structure).

**Interpretation rule (fixed now):** an enclosure failure at stage (b) after a
clean stage (a) is evidence about the bound, and even then I escalate to Main
BEFORE recording any conclusion that a published rung fails to enclose (Main
directive 4; also repo escalation rule). An enclosure failure WITHOUT a clean
stage (a) is evidence about my transcription and must be reported as such, never
as a claim about the bound.

Both rungs get the protocol: rung-1 = Alman25 W1.00_2.371339.mat → ω ≤ 2.371339
(priority: it is the immediate predecessor and the number everyone cites), then
rung-2 = VXXZ24 K100_2.37155181.mat → ω ≤ 2.371552. DWZ23's MCOS-opaque .mat may
be unreadable without MATLAB; if so I say so and stop there rather than guessing

## Addendum A3 (2026-08-30T01:35Z, per Main directive — the feasibility gap)

Stage (a) measured, at the exact released point of rung 1 (Alman25
W1.00_2.371339.mat): max inequality violation c_max = 1.137e-10, max equality
violation ceq_max = 2.390e-11, Schönage-line quantity
(retained + M·ω − 2^(ℓ−1)·ln(q+2), with ℓ=3, q=5) = −6.2e-15 measured in
float64. The Schönage line is a **`c`-type constraint written as
`target − value ≤ 0`, so `value − target ≥ 0` is required; the measured −6.2e-15
is a VIOLATION of magnitude 6.2e-15**, i.e. the released float point sits on
the wrong side by half-ulp float noise. All three residuals are small but
NONZERO: the released point is NOT exactly feasible, and a 1e-9 refine
tolerance is a numerical convergence criterion, not evidence of an exactly
feasible point. Treating these residuals as zero would certify nothing.

RESOLUTION CHOSEN (re-scoping gate B accordingly): **(a) SLACK ABSORPTION**.
The re-scoped gate-B deliverable is a certified statement of the form

    ω ≤ 2.371339 + ε_subject_to_feasibility,

where each retained-count R and the matrix size M are lowered by the certified
outward-rounded cost of absorbing the constraint slack (each c-constraint
`g(θ) ≤ 0` violated by δ is absorbed by lowering the corresponding R by the
dual-priced least-amount, and the max-entropy penalties H^max − H are charged
their Lemma-1 ε already), and ε is EXPLICIT and outward-rounded. If ε comes out
≈ 1e-8 or smaller, this is the first rigorous enclosure of the rung (up to the
absorbed slack), and it is reported with ε never hidden. Exact repair (option
b) is attempted only if slack absorption gives ε > 1e-6; honest negative
(option c) is the fallback with the explicit least-ε statement.

The stage-(b) interval evaluator charges, per hashing constraint, the
interval-enclosed quantity `num_block_lower − P_true_upper` (P_true from the
Lemma-1 dual certificates with rigorous ε), and the Schönage line is checked
with retained counts at their certified-feasible lower endpoints and M at its
certified-feasible maximum. The output is thus exactly the "omega ≤ published +
ε" statement of option (a), not an unconditional certified endpoint.

Report format commitment: every gate-B verdict states ω_cert(PASS/FAIL),
the explicit ε, the per-constraint budget of absorbed slack, and whether the
slack was dominated by float noise or by Lemma-1 entropy gaps.

## Gate C — scope (only if B leaves time)

Interval-certified branch-and-bound in a box around the best certified feasible
point from B. Explicit box = conv(v) from gate B. Budget: same 10⁶-leaf cap per
subdivision, single process, `nice -n 10`. Output: either a certified local-optimality
cap (all neighbours in the searched box ≥ certified bound of the found point) or a
certified improvement certificate. SEPARATE evidence track: MACHINE-VERIFIED only
if a complete VeriPB/branch-bound proof is materialized; otherwise
COMPUTATIONAL-EVIDENCE.

## Anti-hindsight clause

The thresholds, precision, budget, certificate format in this file were fixed
BEFORE any gate-B numeric run. Gate A finding fixed before any from-scratch
optimization attempt of ours (none planned for gate A). No change to the above is
permitted without an addendum below recording the change and its reason.

## Addenda

- (none yet)
