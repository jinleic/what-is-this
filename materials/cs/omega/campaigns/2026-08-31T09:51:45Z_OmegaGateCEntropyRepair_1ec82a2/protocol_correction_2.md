# Protocol correction — Main five-point audit response (2026-08-31)

Status: STATIC ONLY. Zero entropy solves and zero objective evaluations before
or after this correction.

## Implementation of the five held items

1. **Frozen dependency loading (item 1).** `gate_c_entropy_repair.py` now
   resolves `gate_c_candidate_witness_dependency`, `interval_core`, and
   `vxxz24_float` from this campaign's checksummed `.asrun` copies via
   explicit `importlib` loading (`SourceFileLoader`; no live
   `omega/src` import). At startup it asserts the SHA-256 of the frozen
   frozen `interval_core.py.asrun` and `vxxz24_float.py.asrun` against the
   corresponding `omega/src` file (the witness dependency has no live
   counterpart), and writes each dependency's frozen SHA-256 into
   `protocol_guards.frozen_checksummed_dependencies` with
   `live_vs_frozen_hash_assertion: true`. The pending run command executes
   this audited source itself; no live-source dependency remains.

2. **Exact-rational primal / full-row feasibility (item 2).** New helpers:
   - `feasible_primal_lower`: exact Fractions test on ALL B rows
     (`sum(x)=1`, every x >= 0, every `B x` exactly equal, and the stored
     `dist_max` probabilities compared by exact Fraction equality). The
     replaced `b_full` is returned so the caller can still use the interval
     dual; failure falls back to H(p).
   - `exact_rational_primal`: for an independent row set, solve the
     underdetermined exact system via SymPy `gauss_jordan_solve` with integer
     (0/1) backfill of the remaining free symbols, reject anything
     non-smooth, then check the exact row equations exactly. This is the
     exact near-dual primal projection; all rejections are recorded in
     `entropy_certificates` metadata as an honest FALLBACK rather than
     asserted feasibility.
   - Infeasible-overlap control: a per-branch plant constructs a max-vector
     with an explicit `[-1,0]` interval on every coordinate and asserts it is
     REJECTED (`status != "stored dist_max exact-feasible"`); this directly
     addresses the old interval-inclusion bug (overlap accepted what the
     exact test rejects).

3. **Entropy-only suppression of Lagrange equalities (item 3).**
   `wrap_class_method` monkey-patches, on the frozen checksummed module's
   *live class objects* (the ones `Workspace` actually instantiates after
   loading), both `Part.lagrange_constraints` and
   `GlobalStage.lagrange_constraints` to return `[]` once per process
   (`_repair_wrapped` guard). The two structural rows at
   `Workspace.evaluate` (line 1190, `ceq_viol.append(self.ret_comp[(r,2)])`
   for `r = 1..2`) are NOT touched, so `equality_abs` in the defect term
   still charges exactly those two structural equalities.

4. **Repaired penalties inside evaluate_post for constraint alignment
   (item 4).** `Part.evaluate_post` and `GlobalStage.evaluate_post` are
   wrapped so that after the frozen methods run and populate
   `num_block_contribution` / `num_block` normally, the wrapper overwrites
   `hash_penalty_term[r]` with the repaired Arb interval
   `gap * scale` (scale = region_prop x part_frac for Part). The manual R
   loop then reads the SAME object fields via
   `patched_part_penalty`/`patched_gm_penalty`, so `Workspace.c_viol`
   hashing constraints and the manual R aggregation are the identical
   intervals. The wrapper asserts a per-branch key signature
   `(point, kind, region, shape)` for uniqueness, so no branch can be
   skipped or double-charged.

5. **KKT-only mandatory plant (item 5).** The counterfactual is now purely a
   validity check: perturb the glyph-multiplier vector by the exact dyadic
   `+1/8` on its first in-row slot, assert
   `planted_kkt_upper > 2^-32` (hard requirement), and log the expected
   fragmentation — the planted `planted_upper.hi > unplanted upper.hi` is
   recorded as `strictly_worse_diagnostic` but is NOT part of `caught` (it
   is optimizer-dependent and would let a certified-validity failure slip
   through if over-tightened).

## Audited-source replacements in place

- `omega/src/gate_c_entropy_repair.py` — the audited source, committed as a
  static refactor snapshot. The pending run command targets this file.
- New `protocol_static_v3.json` — refreshed exact-guard record from the
  prepare-only smoke run (`kkt_limit`, per-point coordinate exact
  coordinates and frozen checksums recorded).
- `checksums.sha256` — refreshed below after adding this file + the v2/v3
  static outputs; the audited source is also executed from the campaign
  path via the `gate_c_entropy_repair.py` committed identity in checksums.

## Not yet done, unchanged

- No entropy solve and no objective evaluation: the pending run command is
  still gated on the shared heavy slot.
- `equalities` handling in the frozen `interval_core.asrun`'s
  `maximum_absolute` remains the single consumed path for the two
  structural equality rows; no change was made there.
