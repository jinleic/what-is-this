# Cycle preregistration — physics/qlops — 20260908T193320Z_ba7b32

hmz FREEZE artifact and campaign prereg for this cycle's single bounded
campaign. Written 2026-09-08T20:1xZ, AFTER the bounded F-Z4 mechanism
diagnosis (sub-gate 1, no author-count sampling) and BEFORE
`campaign.py init` and before any Revision-8 sampling.

## Gate question

Owner directive: resolve finding F-Z4 WITHOUT depending on the authors, via
the Revision-8 route. F-Z4 (frozen run
`20260908T174922Z_dc3b41ff_d04f0ae58ca0`, REJECTED): shipped grown
`832_text_{0.0008,0.0006,0.0004,0.0001}.stim` differ from the pinned-driver
rebuild by an identical topology delta. The question: (i) what is the
mechanism of the split — (a) serialization/flattening artifact, (b)
driver-version drift, or (c) genuinely different circuits — decided from
structure alone; and (ii) under the preregistered Revision-8 amendment
(rebuilt-oracle equality for the four split cells: the pinned-driver
rebuild IS the reference; shipped artifacts demoted to provenance), does
the full 12-cell / 71,000,000-shot author schedule reproduce under the
pinned driver?

## Sub-gate 1 result (completed BEFORE this prereg's run half)

Mechanism verdict: **(b) driver-version drift.** With the pinned driver
(commit `1b59e22`, adapter exec path, on mini-0): all four split cells show
an identical bounded delta — same 968-qubit lattice and coordinate sets;
after `QUBIT_COORDS` relabeling, all 841 DETECTOR definitions and all 3
OBSERVABLE_INCLUDEs byte-identical; the delta is confined to 5 of 36
tick-rounds (rounds 10–14) where three ancilla qubits are
reset/wired/measured on a different schedule (net +13 CX, −3 R, −2 X_ERROR
shipped-side, DEPOLARIZE2 placement matching); DEMs share 14,307/14,466
lines (Jaccard 0.989); interaction graphs WL-equal; cross-p
noise-normalized comparison: six rebuilds are one generation, the four
split shipped files another. Eliminated: (a) — instruction multisets,
measurement streams, and DEMs genuinely differ; coordinate relabeling does
not reconcile them. Eliminated: (c) — lattice/detectors/observables
identical, delta bounded and systematic. Because DEMs differ, the split
shipped artifacts cannot oracle the pinned driver; because
lattice/detectors/observables are identical, the pinned-driver rebuild is
an honest-broker reference. Evidence: `scratch/diag-report.json`,
`scratch/diag2-report.json`, `scratch/diff_grown_*.txt`,
`scratch/coordrewrite_diff_grown_*.txt`.

## Init-command deviation (recorded before init, not after)

Brief step 3's literal gate slug
`resolve-f-z4-obtain-the-authors-construction-rev` is mechanically
incompatible with the control plane, for the same three reasons documented
in last cycle's prereg: the adapter recomputes the run-id hash12 from its
GATE constant, binds `physics/qlops/pre_statement.md` bytes, and pins the
latest revision heading. Revision 8 (just appended to
`pre_statement.md`) names the replacement gate
`zero-level-author-repro-r8`, adapter version `r8`. The cycle therefore
inits with:

```sh
python3 ../../scripts/campaign.py init \
  --gate zero-level-author-repro-r8 \
  --prereg pre_statement.md
```

This file remains the hmz freeze/cycle-prereg artifact; the machine-bound
prereg is `pre_statement.md` Revision 8. Nothing in this deviation widens
scope.

## Falsifiable acceptance criteria (Revision 8, frozen above)

- Full mode executes the eight-gate smoke battery FIRST; any smoke-gate
  failure ⇒ REJECTED refusal before a single author-count shot.
- 12 cells fresh under Revision 8: 8 shipped-oracle cells (all six ungrown
  + grown@{0.001,0.0002}; shipped-vs-rebuilt flattened inequality ⇒
  REJECTED refusal before sampling that cell) and the 4 F-Z4 split cells
  under rebuilt-oracle semantics (pinned-driver rebuild is the sampled
  reference; shipped artifact provenance-only; factual inequality
  recorded, never refused).
- 24 comparisons (12 LER + 12 acceptance) against the author rows embedded
  in the pinned sources, combined analytic binomial sigmas, familywise 1%
  over 24 dependent comparisons, Bonferroni-derived |z| ≤ 3.53 gate.
- Preregistered verdict mapping (adapter `verdict_from_zs`, unchanged):
  any |z| ≥ 5 or ≥ 3 values beyond 3.53 ⇒ FROZEN-NEGATIVE; otherwise
  undefined statistics or 1–2 mild failures ⇒ FROZEN-INCONCLUSIVE; all 24
  defined |z| ≤ 3.53 ⇒ FROZEN-CERTIFIED.

## Verdict mapping for `campaign.py close --verdict`

- Adapter completes full mode, candidate FROZEN-CERTIFIED (12/12 under
  Revision 8) ⇒ close **FROZEN-CERTIFIED** — the reproduction completes
  with the amendment documented (certifies the executable d3/d9 artifact
  ONLY, with the four shipped split artifacts documented as F-Z4
  provenance).
- Adapter completes full mode, candidate FROZEN-NEGATIVE ⇒ close
  **FROZEN-NEGATIVE** — the split is real at the statistics level; the
  mechanism ((b) driver-version drift) is named.
- Adapter completes full mode, candidate FROZEN-INCONCLUSIVE ⇒ close
  **FROZEN-INCONCLUSIVE**.
- Adapter refuses (exit 1, REJECTED) before sampling ⇒ close **REJECTED**.
- Run interrupted before `results/analysis.json` exists (cycle wall-clock
  ceiling 3.5 h reached) ⇒ stop at the last completed per-point artifact,
  freeze the partial state, close **FROZEN-INCONCLUSIVE**, and report the
  exact completed fraction (cells and shots) — partial progress is a
  result, never a certification.
- Adapter/process crash (no refusal, no analysis) ⇒ close **CRASHED**.

## Scope boundary

- IN: the F-Z4 mechanism diagnosis (done, sub-gate 1); the Revision-8
  amendment of `pre_statement.md`; the adapter r8 amendment (oracle
  semantics for the four split cells; gate/version/revision-marker bump;
  test-file contract updates); ONE fresh 71m-shot/12-cell full run under
  Revision 8; its audit and close.
- OUT: paper d=7 and physical c≈300 (remain NOT-REPRODUCED regardless of
  outcome; no d=7 reconstruction); Gate A/B conclusions (unchanged);
  contacting the authors or fetching any other repo revision (the owner
  directive forbids depending on the authors); any modification of the
  R5/R6/R7/REJECTED frozen run dirs; any change to the source tar;
  reinterpretation of the frozen REJECTED run (its 7 completed cells are
  NOT adopted; the Revision-8 schedule is run fresh end-to-end);
  physics/qldpc-dec and other targets.
- Compute placement: ALL remote compute on `mini-0.local` ONLY (shared
  with the GB9 decoder job: thread caps ≤ 6, `nice -n 10`, launched
  detached via `nohup caffeinate -s … &`). Control plane strictly local
  (campaign.py, state.json, this workspace). Inputs copied out; results
  copied back into the run dir before freeze.

## Exact verification command

On `mini-0.local`, from the replicated target tree
`~/qlops-r8/qlops` (r8 adapter + test file + Revision-8 prereg + the
pinned tar, all sha256-verified against the local copies before launch):

```sh
OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 \
NUMEXPR_NUM_THREADS=6 VECLIB_MAXIMUM_THREADS=6 \
/usr/bin/nice -n 10 caffeinate -s ~/qlops-r8/.venv/bin/python \
  src/zero_level_repo_repro.py \
  --source-tar source_tar.gz \
  --run-dir campaigns/<printed-run-id> \
  --mode full
```

Audit re-verification (fresh, post-run): `sha256sum -c sha256s.txt` inside
the frozen run dir; source tar sha256 must remain
`d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338`;
bitwise replay of one rebuilt-oracle cell (grown@0.0008, seed 2202)
against its stored row.
