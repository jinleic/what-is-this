# PRE-STATEMENT — OmegaGateCEntropyRepair v11 (Fable-audit repair campaign)

Written 2026-08-31T11:30:16Z (D2-D5 + perimeter addendum) before any v11
evaluation. PARENT: `..._v10_4cdd3c85`; GRANDPARENT: v9 (owner-accepted but
real-run ABORTED), preserved untouched.

## Acceptance revocation (external Fable audit)

v10 acceptance was REVOKED by Main citing external Fable-audit production
defects beyond the (correct) interface math. v9 and v10 are preserved
untouched as rejected/aborted history. This campaign stages v11 and holds.

## D2 exact moved dyadics (main defect of v10)

`exact_dyadic_fraction` used `Fraction.from_float(float(mid))`, which is
exact for the ROUNDED binary64 only and silently changed >53-bit exact
Arb dyadics. v11 replaces it with exact finite Arb point extraction:
`arb_exact_fraction` pulls mantissa/exponent from the ball via `man_exp()`,
requires integer mantissa and zero radius, then reconstructs and requires
byte-exact endpoint equality.  Consequences (MACHINE-VERIFIED in the v11
mock):
- a 71-bit dyadic `(2^70+1)/2^70` now round-trips PASSES (old converter
  silently produced 1; recorded as `over53bit_changed_by_old_float`);
- non-point intervals still hard abort (`d2_nonpoint_hard_abort_kept`).

A prepare-only guard now verifies ALL 45 coordinates at EVERY fixed rung
by comparing the exact `center +/- step` Fractions (the exact construction
matching W.get) against the preregistered `coordinates` in
`protocol_guards.steps_and_exact_coordinates` via `Fraction ==` equality.
This structural guard is the only point-consuming preparation and consumes
no objective evaluation.

## D3 interval enclosure fixes (canonical live + v11 frozen + correction)

- `from_two` upper endpoint fixed from `b.lower()` to `b.upper()`.  The
  previous b.lower() truncated the top of a nonzero-radius ball (observed:
  old upper after ln(3/4) is 1 ulp SHORT of `b.upper()`). Regression plant:
  a nonzero-radius Arb ball must pass upper containment — the OLD code
  fails, the NEW code passes (`d3_from_two_upper_containment`).
- `ln_iv` nonpositive-fallback changed from `arb(-10)**40` (= NEGATIVE
  1e40, which nonsensically claims log could be a huge negative) to
  `-(arb(10)**40)` (a very large NEGATIVE magnitude, correct semantics).
  Plant: value is strictly negative (`d3_ln_iv_fallback_negative`).
- 1/e computed via Arb `exp` (not binary64 `np.e`); containment plant
  asserts `1/e` is enclosed by the Arb-derived bounds
  (`d3_one_over_e_arb_containment`).
- `if not __debug__: raise RuntimeError` guard inserted so asserts cannot
  silently strip under `python -O` in production.

## D4 frozen perimeter (no live-source execution)

- `alman25_float.py.asrun_v11` is a frozen checksummed copy of the live
  `alman25_float.py`; it is loaded into `sys.modules` BEFORE
  `interval_core` (which imports alman25_float), preserving the
  IntervalTree default path through frozen code only.  No deletions.
- v11 `interval_core` `sys.path.insert` now points to its own campaign
  directory (via `dirname(abspath(__file__))`), not the live `omega/src`.
- The frozen `gate_c_candidate_witness_dependency.py.asrun_v11` has had
  its `SRC`/`sys.path.insert` rewritten to `os.path.dirname(os.path.abspath(
  __file__))` (its own campaign dir) and imports now resolve the already
  hash-verified `sys.modules` entries.
- Artifact grep confirms ZERO `/omega/src` hits in any executable frozen
  dependency. CANONICAL LIVE FILES (cs/omega/src/interval_core.py and
  cs/omega/src/alman25_float.py) keep their NORMAL sys.path/pth values
  UNCHANGED — v11 freezes specific frozen copies; live paths untouched.

## D5 inputs (local copies + hard-coded hash verification)

Copies into v11 campaign and verified by SHA-256 BEFORE reading:
- `K100_2.37155181.mat`
  (verified: f15684a43d0193fdd2c1ceba9e1001192f7b2bad36f0ca715636dcecf9e66eba)
- `kernel_basis_V_exact_rational.json`
  (verified: 7d69c1dd4c9355dfe6d6b3a762658ba10d93e4bcaea3b240ebbbd811de3eba60)
- `margin_matrix_M.npy`
  (verified: 423dbafd7c530ec7256ddf68cafdc15e143503ed9358c39bb81e4bc66b074b4e)

`W.MAT = V11_MAT` and `W.CAMP0 = HERE` (v11-local) are set BEFORE any
prepare/evaluate.  Any digest mismatch hard-aborts rather than loading.

## Cache-risk fix (CERT_CACHE dual key)

The dual certificate cache key now includes the residual system:
`(residual_matrix.shape, residual_matrix.tobytes(), tuple((lo, hi) for
residual_b))` in addition to the previous key components, so cached KKT
cannot cross residual systems.

## Frozen command (on release)

    cd cs && nice -n 10 env OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
      ./.venv/bin/python \
      omega/campaigns/2026-08-31T11:30:16Z_OmegaGateCEntropyRepair_v11_d2d3d4d5/gate_c_entropy_repair.py.asrun_v11 \
      omega/campaigns/2026-08-31T11:30:16Z_OmegaGateCEntropyRepair_v11_d2d3d4d5/results_v11.json \
      > omega/campaigns/2026-08-31T11:30:16Z_OmegaGateCEntropyRepair_v11_d2d3d4d5/run_v11_stdout.log 2>&1

Five threads=1, nice 10, alarm 900s, fixed point order base →
R/{16,8,4,2,1}, first-strict-stop.  `results_v11.json` + additive post-run
manifest; no pre-run history rewrite; no deletions (pycache preserved,
disclosed as non-evidence).

## ZERO real point compute on this campaign

- `py_compile`
- `--prepare-only` (constructs the real Workspace, runs the structural
  real-object guard, and stops — NO evaluate_point, NO entropy solve,
  NO objective evaluation)
- `--mock-self-test` (mock objects only; the D2/D3 plants are in this
  self-test, all static; no real points)

## Rule 7 scope sentence

Unchanged from parent: swept set remains base plus five fixed steps
R/{16,8,4,2,1} on one exact integer direction d, only the 45 region-0
glob dist[0] coordinates of the frozen VXXZ24 `K100_2.37155181` vector
under the transcribed single-p_comp max_level=3 q=5 program, preserving
mass via sum(new)-sum(base)=0, generalized-entropy/absorbed-defect
objective at observed mass m (never normalized to 1).  Nothing else
searched.  Even a strict result establishes only STRICT DESCENT on the
named fixed-m affine ray at the released binary64 center.  No paper or
record claim.

## Owner correction and final static refreeze

This append-only section supersedes inconsistent staging text above.  The
staging agent's reported v11.3 still converted each moved coordinate through
binary64 and therefore did not exercise the actual `W.get` representation.
The owner replaced that guard with the exact production construction under
300-bit Arb: `IC.const(center) +/- IC.const(step)`, followed by
`exact_dyadic_fraction` and exact equality against the preregistered Fraction.
The final prepare-only record passes **270/270** comparisons: 45 base
coordinates plus 45 coordinates at each of five rungs.  The guard record key is
`protocol_guards.d2_all_rung_converter_guard`; the older
`steps_and_exact_coordinates` wording above is not a JSON path.

The old `ln_iv` expression `arb(-10)**40` is positive because the even power
applies to `-10`; the corrected `-(arb(10)**40)` is negative.  The final
interval core also hulls all four multiplication/division candidate endpoints,
uses outward endpoint pairs for non-dyadic rational entropy inputs, and tests
the entropy maximum with pure-Arb `exp(-1)`.  Owner regression
`owner_interval_regression_v11.json` passes the old-fails/new-passes
`from_two` case, multiplication and division hulls, the $1/e$ entropy maximum,
Fraction $1/3$, the negative-log fallback, and the production `__debug__`
gate.

Canonical live `cs/omega/src/interval_core.py` and
`cs/omega/src/gate_c_entropy_repair.py` were changed and are byte-identical to
the final release files.  Existing campaign aliases were preserved and
synchronized byte-identically; nothing was deleted.  The final release source
hash is
`ca6be84a18cf5e63eeec6fd664ff9e4dcd3e70de90966441996ce314524a7581`;
the interval-core hash is
`ce70f95922f3d5098e892e1b7047570ebb5267826b51f2776c91b56d9934adb5`.

The command block above contains a malformed combined environment assignment.
The authoritative frozen command is the `manifest.json` command, with five
separate assignments:

    cd cs && nice -n 10 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
      ./.venv/bin/python \
      omega/campaigns/2026-08-31T11:30:16Z_OmegaGateCEntropyRepair_v11_d2d3d4d5/gate_c_entropy_repair.py.asrun_v11 \
      omega/campaigns/2026-08-31T11:30:16Z_OmegaGateCEntropyRepair_v11_d2d3d4d5/results_v11.json \
      > omega/campaigns/2026-08-31T11:30:16Z_OmegaGateCEntropyRepair_v11_d2d3d4d5/run_v11_stdout.log 2>&1

Final static evidence is additive:
`prepare_only_v11_owner.json`, `mock_self_test_v11_owner.json`, and
`owner_interval_regression_v11.json`.  No objective point, entropy solve, or
production sign was evaluated.

`results_v11.json` was already populated by an agent prepare-only invocation
from a rejected static iteration.  It is preserved as stale static evidence and
must not be overwritten.  The authoritative manifest therefore reserves
`results_v11_production.json` and `run_v11_production_stdout.log` for the first
production release; these names supersede both command blocks above.

## Owner output-perimeter correction

The `ca6be84a…` release source above still accepted a default/arbitrary output
path and replaced pre-existing output or temporary files. It is revoked before
production. The mathematical certificate path is unchanged, but the accepted
release source is now
`cb86fe9e0a82b32be00f82bed6330c211e719936938ca00e5af79c60518b625b`.
It requires an explicit output directly in this frozen campaign, rejects a
pre-existing output or `.tmp` path at startup, creates each checkpoint temp
exclusively, and fsyncs the campaign directory after atomic rename. The frozen
alias and canonical live source are byte-identical to these release bytes.

The authoritative manifest command also enables shell `noclobber` before
redirecting `run_v11_production_stdout.log`. The older static outputs remain
historical evidence for the unchanged mathematical core, not verification of
the corrected output perimeter. Fresh `cb86fe9e…` prepare-only and mock checks
are required before production.

Fresh static checks of `cb86fe9e…` pass without a real point evaluation.
`prepare_only_v11_owner_perimeter.json` repeats the 270/270 exact-coordinate
guard; `mock_self_test_v11_owner_perimeter.json` repeats the wrapper, interval,
projection, generalized-mass, and control plants. A launch against preserved
`results_v11.json` refused the pre-existing path and left its SHA-256
`3dfb18d6…` unchanged; an off-campaign `/tmp` output was also refused before
prepare. No `.tmp` artifact remains.

Focused read-only Fable audit of release source `cb86fe9e…`, interval core
`ce70f959…`, and the authoritative manifest returned:
**“PASS — no release blocker found.”** The reviewed claim scope was explicitly
the fixed-m generalized-entropy/absorbed-defect ray evidence only, not a
feasible sum-one construction and not an omega bound. The audit made no edits
and executed no campaign code.

The final released shell command sets `PYTHONDONTWRITEBYTECODE=1`; this
preserves historical `__pycache__` bytes and changes neither the audited
source nor the fixed-m evidence scope.

## Production-v11 abort and bounded raw-anchor diagnostic

The first real v11 release exited `1` after 4.42 seconds at
`base["raw"]` containing `W.FROZEN_RAW`. It preserved
`results_v11_production.json` (`c74cb606…`) and combined stdout/stderr
(`2bd2031b…`). The checkpoint contains only `protocol_guards` and the
`bad_multiplier_plant`; no base checkpoint, rung, or decision exists.

Source reread shows a contract mismatch to test, not a numerical outcome:
`FROZEN_RAW` is the dependency's pre-repair raw-objective anchor, while v11
installs repaired entropy penalties before `evaluate_point`; those penalties
change `total_R` and hence `(target-R)/M`. A bounded diagnostic is registered
to run the same base once from the exact failed source, save to
`raw_anchor_diagnostic_v11.json`, and print the repaired raw interval, old
anchor, repaired `R`, `M`, defect, and full interval to
`raw_anchor_diagnostic_v11_stdout.log`. It runs no rung and makes no
acceptance decision. The failed production artifacts will not be overwritten.

### Raw-anchor diagnostic result and registered rejection diagnostic

The bounded base completed in 4.01 seconds. The repaired raw interval is
`[2.3715518061851393, 2.3715522159185425]`; the old frozen anchor is
`[2.37155383583508, 2.3715538358350807]`, strictly above and disjoint.
The stale containment assertion is therefore inapplicable to the repaired
objective.

The same run also reported
`all_entropy_multipliers_accepted = false`. Removing only the stale anchor
would therefore yield no accepted strict rung. Before freezing a successor,
a second bounded one-base diagnostic is registered under
`entropy_rejection_diagnostic_v11.json` and
`entropy_rejection_diagnostic_v11_stdout.log`. It will retain every rejected
certificate's kind, region, acceptance bit, binary64 multiplier vector, dual
and entropy-gap intervals, and KKT upper bound, plus aggregate rejection
counts. It runs no rung and makes no acceptance decision.

### Entropy-rejection diagnostic result and solver experiment

The base has 192 certificates: 53 accepted and 139 rejected by the exact
`1/2^32` KKT gate. Of the rejected certificates, 136 report optimizer success;
three global certificates hit the 2,000-iteration limit. Rejected outward KKT
upper bounds range from `2.381769915623674e-10` (just above the gate) to
`2.938341911665688e-06`. All 139 rejected records retain exact-mass feasible
primal projections; this isolates the immediate acceptance failure to
multiplier accuracy rather than the mass-feasibility repair.

A bounded one-base solver experiment is registered under
`least_squares_solver_diagnostic_v11.json` and
`least_squares_solver_diagnostic_v11_stdout.log`. For every certificate it
will run the existing L-BFGS-B solver, then use that point to initialize
SciPy `least_squares` on
`B (m softmax(B^T y)) - b`, with analytic Jacobian
`B (diag(q) - q q^T/m) B^T`, tolerances `1e-15`, and at most 4,000 function
evaluations. It records aggregate certified acceptance and solver-path
statistics. It runs no rung and makes no decision.

#### Solver-instrumentation abort and retry

The first solver experiment exited `1` after its one base. `runpy` returned a
namespace distinct from the functions' live `__globals__`; assignment through
the returned mapping therefore made zero replacement calls, and the summary
failed on an empty solver-path list. No result artifact and no solver-math
evidence were created. The failed stdout/stderr is retained.

A corrected retry is registered under
`least_squares_solver_diagnostic_v11_retry1.json` and
`least_squares_solver_diagnostic_v11_retry1_stdout.log`. It installs the same
pre-stated replacement through `entropy_certificate.__globals__` and asserts
the live global identity before `prepare()`. Algorithm, one-base scope,
4,000-evaluation cap, no-rung restriction, and no-decision restriction are
unchanged.

#### Corrected solver diagnostic result

The corrected one-base retry completed in 4.10 seconds and invoked the
replacement exactly 192 times. All 192 certificates passed; all 192
least-squares solves reported success. The maximum binary64 residual was
`7.228612923928862e-11`; the maximum outward certified KKT upper bound was
`1.4456822252392066e-10`, below the unchanged exact gate
`1/4294967296`.

This proves only that the two-stage multiplier solver clears the repaired
base's KKT gates under the unchanged interval checker. It is not rung evidence
and not a decision. A successor must integrate the solver, remove the
inapplicable pre-repair raw anchor, refreeze independently, and rerun every
specified rung.
