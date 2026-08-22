# Research log

All recorded work occurred on **2026-08-11**.  The ordering follows the numbered research sequence; several computations ran concurrently, so this is a logical chronology rather than a fabricated wall-clock serialisation.  Every cycle records the question, method, artifact, outcome, and lesson present in the repository.

## L00 — Source and claim triage

- **Question:** Which historical, numerical, and claimed-solution statements are precise enough to test?
- **Method:** Local source-manifest/full-text review; extraction of formulas, conventions, citations, and invariant predictions.
- **Artifacts:** `notes/literature_matrix.md`, `notes/lit_claimed_solutions.md`, `sources/manifest.yaml`, `results/reference_data/claimed_solutions.json`.
- **Outcome:** Classical exact 2D work, rigorous inequalities, external numerical benchmarks, and three claimed 3D solutions were separated by evidence status.
- **Learned:** Publication and formula complexity are not evidence of exactness; a candidate must predict something not used to choose it.

## E01 — Exact finite-volume foundation

- **Question:** Do brute-force enumeration and transfer propagation implement the same lattice conventions?
- **Method:** Compare complete exact broken-bond polynomials on periodic, open, and mixed-boundary lattices.
- **Artifacts:** `experiments/e01_validate_enumeration.py`, `tests/test_tm_vs_enumeration.py`, `src/ising/exact_enumeration/__init__.py`, `src/ising/transfer_matrix/__init__.py`.
- **Outcome:** 17 coefficient-by-coefficient cross-checks pass.
- **Learned:** Parallel length-2 bonds and dropped length-1 self-bonds must be fixed before any higher-level result.

## E02 — Onsager/Kaufman control

- **Question:** Can the repository reproduce a genuinely solved model at finite volume, thermodynamic limit, and operator-algebra level?
- **Method:** Exact bivariate transfer polynomials, Kaufman's four sectors at 60 digits, Onsager's double integral, Jordan--Wigner identities, matrix logarithms, and generator reconstruction.
- **Artifacts:** `experiments/e02_onsager_2d.py`, `notes/onsager_derivation.md`, `results/onsager/onsager_2d.json`, `tests/test_onsager.py`.
- **Outcome:** Eight finite tori agree to maximum relative error `3.63867276317e-60`; parity-resolved `log V` has maximum off-Majorana-bilinear coefficient `1.246e-15`.
- **Learned:** Boundary sectors and parity resolution are load-bearing; nearby sign and zero-mode conventions fail exact finite integers.

## E03 — First exact finite-lattice series

- **Question:** Can simple-cubic HT/LT free-energy coefficients be derived internally rather than copied from published tables?
- **Method:** Finite-lattice Möbius inversion with exact `Fraction` arithmetic and exact transfer polynomials; 2D Onsager control.
- **Artifacts:** `experiments/e03_flm_series.py`, `notes/flm_derivation.md`, `results/series/sc_ht_free_energy.json`, `results/series/sc_lt_free_energy.json`, `tests/test_series.py`.
- **Outcome:** Initial simple-cubic series reached HT `v^16` and LT `x^20`; the square-lattice control matched Onsager through `v^24`.
- **Learned:** Connected-cluster box budgets and exact boundary normalisations, not numerical fitting, determine the available order.

## E04 — Lee--Yang finite zeros

- **Question:** Do exact finite 3D field polynomials obey the Lee--Yang circle theorem, and can the available sizes determine an edge exponent?
- **Method:** Exact joint DOS/field polynomials, palindrome and zero-field identities, 100-digit roots, 2D controls, and finite-size edge fits.
- **Artifacts:** `experiments/e04_lee_yang.py`, `notes/lee_yang.md`, `results/lee_yang/lee_yang_analysis.json`, `tests/test_lee_yang.py`.
- **Outcome:** Unit-circle/root-count checks pass on 13 lattices through `4x4x4`; the largest stored 64-site radial residual is `9.1044993726259e-92`.  The edge-exponent estimate is rejected because the same procedure is unreliable on 2D controls.
- **Learned:** Exact finite theorem verification is strong; a scaling exponent from a few small boxes is not.

## E05 — Rigorous critical-coupling interval

- **Question:** Can `K_c` be enclosed without using the numerical benchmark as an input?
- **Method:** Exact self-avoiding-walk counts and correlation domination for the lower bound; infrared domination and the simple-cubic Watson integral for the upper bound; directed rounding.
- **Artifacts:** `experiments/e05_kc_bounds.py`, `proofs/kc_bounds.md`, `results/bounds/kc_bounds.json`, `tests/test_bounds.py`.
- **Outcome:** `K_c` is certified in `[0.2074277114992039908436804465100887760027, 0.2527310098586630030260020266135701299926]`.
- **Learned:** The problem-stated Watson gamma formula with denominator `4*pi^3` is too large by 8; the correct denominator is `32*pi^3`.

## E06 — Falsification of claimed exact solutions

- **Question:** Do the Zhang 2007 and Degang Zhang 2021 proposals pass invariant critical, symmetry, and series checks?
- **Method:** High-precision critical roots, exact LT magnetisation coefficients, anisotropic permutation tests, exact HT expansion, and finite-torus winding projection.
- **Artifacts:** `experiments/e06_falsify_claims.py`, `notes/falsification.md`, `results/falsification/falsification.json`, `tests/test_falsification.py`.
- **Outcome:** Zhang's `K_c=0.2406059125...` differs from `0.221654626(5)` by about 3.79 million quoted standard deviations and its `u^6` coefficient is `-18` instead of exact `+14`.  Degang Zhang's critical line is rotation-dependent and its free energy first fails at `K^4` (`-3` instead of `+3`).
- **Learned:** Low-order graph coefficients and lattice symmetry are cheaper and more decisive than fitting a benchmark.

## E07 — Dynamical Lie algebra of individual terms

- **Question:** Does moving from a 1D layer to a 2D layer destroy the small free-fermion algebra generated by individual `X_i` and `Z_iZ_j` terms?
- **Method:** Exact modular Pauli-Lie closure for chains, rings, and open grids.
- **Artifacts:** `experiments/e07_dla.py`, `results/dla_table.json`, `results/e07_dla.log`, `proofs/algebraic_obstruction.md`.
- **Outcome:** Open chains follow `n(2n-1)`, rings `2n(2n-1)` in the tested range; grid dimensions are `56`, `1056`, `16256`, `65535` for `2x2`, `2x3`, `2x4`, `3x3`.
- **Learned:** A genuine 2D layer rapidly fills a large even-parity support space, unlike the free-fermion controls.

## E08 — Algebra generated by the two sums

- **Question:** Does the much smaller Onsager algebra generated only by `A=sum X` and `B=sum ZZ` remain small in 3D layers?
- **Method:** Exact Pauli-Lie closure of the two summed generators over two primes.
- **Artifacts:** `experiments/e08_onsager_algebra.py`, `results/onsager_algebra_dims.json`, `results/e08_onsager_algebra.log`, `proofs/algebraic_obstruction.md`.
- **Outcome:** Open-chain dimensions are `n^2`, ring dimensions `3n-1` in the tested range; open grids give `11`, `263`, and `2952` for `2x2`, `2x3`, `2x4`.
- **Learned:** The 3D-layer obstruction is visible even after summing the local terms, but dimension must be measured separately from term-DLA size.

## E09 — Exact Ising/gauge duality

- **Question:** Is the simple-cubic Ising model self-dual, or dual to a distinct gauge theory with nontrivial sectors?
- **Method:** Exact GF(2) cell-complex homology, gauge-orbit DOS, surface sectors, dual-spin enumeration, and 2D four-sector controls.
- **Artifacts:** `experiments/e09_duality.py`, `proofs/duality_3d.md`, `results/duality/duality.json`, `tests/test_duality.py`.
- **Outcome:** 27 exact integer checks pass.  Three-dimensional Ising is dual to a `Z_2` gauge model, not to itself; naive single-sector torus identities fail until all homology/twist sectors are included.
- **Learned:** `sinh(2K)sinh(2K*)=1` does not imply a 3D Ising self-dual critical point.

## E10 — Local conserved-charge search

- **Question:** Does the 2D layer possess a range-growing tower analogous to 1D free-fermion charges?
- **Method:** Exact finite-field kernels for translation-invariant connected-support Pauli ansätze at generic couplings, with 1D controls.
- **Artifacts:** `experiments/e10_conserved_charges.py`, `notes/integrability.md`, `results/integrability/conserved_charges.json`.
- **Outcome:** The 1D control has nontrivial dimensions `1,3,5,7`; the 2D layer has only `span{I,H}` through range 4 in an ansatz of dimension 2461.
- **Learned:** The bounded local Pauli-charge mechanism is absent, while unrestricted full matrix commutants are too generic to diagnose integrability.

## E11 — Tetrahedron-equation diagnostic

- **Question:** Does the isotropic bond-dimensional-2 Ising tensor satisfy the tested constant tetrahedron equation?
- **Method:** Exact residual polynomials and metric-preserving orthogonal/monomial gauge checks.
- **Artifacts:** `experiments/e11_tetrahedron.py`, `notes/integrability.md`, `results/integrability/tetrahedron.json`.
- **Outcome:** A physical 3D specialization has a nonzero exact residual; permitted gauges do not remove it.
- **Learned:** An unrestricted GL(2) change alters the copied-spin bond metric, so it cannot be counted as the same delta-bond tensor model.

## E12 — Simon--Lieb finite-box bound

- **Question:** Can exact finite-box boundary correlations independently improve the rigorous lower endpoint?
- **Method:** Exact rational correlation polynomials, high-precision roots, conservative forward-error and decimal rounding certificates.
- **Artifacts:** `experiments/e12_simon_lieb.py`, `proofs/simon_lieb_bound.md`, `results/bounds/simon_lieb_bounds.json`, `tests/test_simon_lieb.py`.
- **Outcome:** `K_c>0.20517823158099` is rigorous but weaker than the SAW lower endpoint.
- **Learned:** An independently valid method need not improve the best bound; sharpening the criterion matters more than simply quoting a larger box.

## E13 — Minimal graph obstruction and claw correction

- **Question:** Is an induced claw equivalent to explosive growth of the two-sum algebra, and how sharply does one extra bond change a path?
- **Method:** Exact two-prime closures on paths, cycles, stars, complete graphs, ladders, and path-plus-rung graphs.
- **Artifacts:** `experiments/e13_minimal_obstruction.py`, `results/minimal_obstruction.json`, `results/e13_minimal.log`.
- **Outcome:** One extra bond gives dimensions `470`, `3426`, `6397` for `n=6,7,8`.  But clawed `K_(1,3)`, `K_(1,4)`, and `K_4` have only `24`, `33`, and `15` dimensions.
- **Learned:** **“claw iff exponentially large algebra” is false.**  A claw is a local Majorana/DG obstruction, not a universal dimension theorem.

## E14 — First Dolan--Grady computation

- **Question:** Which Dolan--Grady relation fails on a square-grid layer, and what support causes the failure?
- **Method:** Exact phase-free Pauli commutators on chains, rings, stars, grids, and a cubic layer.
- **Artifacts:** `experiments/e14_dolan_grady.py`, `results/dolan_grady.json`, `proofs/dolan_grady_defect.md`.
- **Outcome:** DG1 holds in every tested case; DG2 fails exactly at degree-3-or-higher vertices and contains a quartic component.
- **Learned:** The failure is not numerical and is localized at layer-graph branching.

## E15 — Local structure of the DG residual

- **Question:** How does the defect depend on vertex degree?
- **Method:** Exact star calculations and grouping of Pauli words by local support.
- **Artifacts:** `experiments/e15_dg_local_structure.py`, `results/dg_local_structure.json`, `proofs/dolan_grady_defect.md`.
- **Outcome:** The residual is a degree-linear term plus a fixed sum over neighbor triples; its quartic multiplicity is `sum_v C(deg(v),3)`.
- **Learned:** The exploratory conjecture that degree `d` creates a `(d+1)`-body all-neighbour term is false; body order remains quartic.

## E16 — DG closed-form theorem check

- **Question:** Does the local formula agree term-for-term on arbitrary finite graphs?
- **Method:** Exact integer comparison of every Pauli coefficient on 16 graphs, plus the local commutator proof.
- **Artifacts:** `experiments/e16_dg_theorem_check.py`, `results/dg_theorem_check.json`, `proofs/dolan_grady_defect.md`, `tests/test_gates.py`.
- **Outcome:** The closed form is machine-verified on all 16 cases and proved locally; the quartic defect vanishes exactly for maximum degree at most 2.
- **Learned:** The precise zero-residual quantifier needed a later disconnected-graph audit; see A01.

## E17 — Symmetry-adapted DG deformation search

- **Question:** Can local invariant deformations cancel the quartic defect and restore both DG relations?
- **Method:** Exact orbit bases, rational rank systems, finite candidate floating constants, and selected Gröbner classifications.
- **Artifacts:** `experiments/e17_dg_deformation.py`, `proofs/dg_deformation_nogo.md`, `results/deformation/dg_deformation_search.json`, `tests/test_dg_deformation.py`.
- **Outcome:** An explicit F3 deformation cancels DG2 on the `3x3` torus, but has a DG1 witness: coefficient `3` in `ad_A'^3(B)` where the matching `ad_A'(B)` coefficient is `0`.  Several B-only/edge-weight families are exactly ruled out.
- **Learned:** Quartic cancellation is possible but insufficient for an Onsager algebra; unrestricted simultaneous A+B deformation remains open.

## E18 — Modular-CRT series extension

- **Question:** Can exact series exceed the 62-site int64 transfer ceiling?
- **Method:** Modular transfer over three 31-bit primes, CRT reconstruction, larger canonical boxes, and prefix self-consistency checks.
- **Artifacts:** `experiments/e18_series_extend.py`, `src/ising/transfer_matrix/crt.py`, `results/series/extended_sc_ht_free_energy.json`, `results/series/extended_sc_lt_free_energy.json`, `tests/test_series_extend.py`.
- **Outcome:** HT reaches `v^20`; LT reaches `x^28`; every previous coefficient remains unchanged.  Direct GF(2) connected-cluster enumeration independently agrees through `v^10`.
- **Learned:** CRT removes overflow, but box count and cross-section still set a hard resource boundary.

## E19 — Short-series singularity and D-finiteness analysis

- **Question:** What can ten nonzero HT coefficients say without benchmark tuning?
- **Method:** Untuned Dlog--Padé/differential approximants at `v^16,v^18,v^20`; exact Euler-ODE fits with three held-out coefficients and an `exp(z)` control.
- **Artifacts:** `experiments/e19_series_analysis.py`, `notes/series_analysis.md`, `results/series/extended_series_analysis.json`.
- **Outcome:** `K_c=0.2229380779607 +/- 0.0050376391291` is a numerical stability envelope containing the external benchmark.  The exponent is weak (`alpha=0.34698 +/- 0.34245`).  No tested low-budget ODE survives holdout; the control is found exactly.
- **Learned:** The series locates the singularity better than the exponent and cannot support a credible exact-constant search.

## E20 — Depth-graded algebra growth

- **Question:** Can finite exponential growth be certified without calling a fit a theorem?
- **Method:** Exact point-group orbit ranks and grading blocks over two primes; depth filtration `D_k`; recorded resource walls.
- **Artifacts:** `experiments/e20_algebra_growth.py`, `notes/algebra_growth.md`, `results/algebra_growth/profiles.json`, `results/algebra_growth/exact_dimensions.json`, `tests/test_algebra_growth.py`.
- **Outcome:** Controls reproduce `n^2` and `3n-1` through `n=32`; grid values/lower bounds include `11,263,2952,8034`.  `D_(2L)>=2^L` is certified only for `L=2..6`; `2x5` reaches `D_20>=7183` at the time wall.
- **Learned:** These are rigorous finite modular lower bounds.  No all-length exponential theorem exists yet.

## E21 — General tridiagonal relation decision

- **Question:** Could the undeformed 3D pair satisfy a general tridiagonal/Askey--Wilson relation even though DG fails?
- **Method:** Exact rational linear systems for all affine parameters; cycles as controls.
- **Artifacts:** `experiments/e21_tridiagonal.py`, `results/tridiagonal.json`, `proofs/tridiagonal_nogo.md`.
- **Outcome:** Cycles have solutions; tested 3D torus/open-grid layers do not.
- **Learned:** The obstruction extends beyond the particular Dolan--Grady normalization.

## E22 — Short tridiagonal certificates

- **Question:** Can each no-solution rank result be reduced to a human-checkable Pauli witness?
- **Method:** Extract minimal inconsistent row subsets from exact systems.
- **Artifacts:** `experiments/e22_tridiagonal_certificate.py`, `results/tridiagonal_certificate.json`, `proofs/tridiagonal_nogo.md`.
- **Outcome:** Short one-/four-row exact certificates were found for the initial cases.
- **Learned:** A convenient single-row pattern looked broader than it really was; E23 tested the quantifier.

## E23 — General certificate audit

- **Question:** Does the proposed one-row witness work on all rectangular layer sizes covered by the no-go?
- **Method:** Exact certificate search on `3x3` through `5x5` and selected rectangles.
- **Artifacts:** `experiments/e23_certificate_general.py`, `results/tridiagonal_certificate_general.json`, `proofs/tridiagonal_nogo.md`.
- **Outcome:** The no-go persists, but the universal one-row wording is false.  It works for all `4x4` strings and half the cases when one side is 4; `3x3`, `3x5`, `3x6`, and `5x5` require four-row witnesses.
- **Learned:** Correct the certificate scope rather than weakening an independently valid theorem.

## E24 — Kac--Ward route

- **Question:** Can a local directed-edge determinant generalize the exact 2D Kac--Ward identity to 3D?
- **Method:** Exact `Z[zeta_8]`/CRT reconstruction in 2D; non-backtracking walk equations and Gröbner elimination in a strict cubic-covariant two-weight slice; dense finite-box reproduction.
- **Artifacts:** `experiments/e24_kac_ward.py`, `proofs/kac_ward_3d_obstruction.md`, `results/kac_ward/kac_ward.json`, `tests/test_kac_ward.py`.
- **Outcome:** Four 2D controls match coefficient-by-coefficient.  The natural unsigned half-angle guess fails at `v^6`; the strict two-weight 3D slice is inconsistent at `v^8`.  The full 30-weight family remains unresolved.
- **Learned:** Planar turning phases are load-bearing, but a scoped scalar no-go cannot be promoted to all local weights.

## E25 — CFT/RG necessary-condition filter

- **Question:** Do candidate free energies have the external three-dimensional specific-heat singularity?
- **Method:** Benchmark-blind curvature scaling with exact 2D logarithmic and synthetic power-law controls.
- **Artifacts:** `experiments/e25_cft_constraints.py`, `notes/cft_constraints.md`, `results/cft/cft_constraints.json`.
- **Outcome:** Both evaluable claimed formulas are logarithmic (`alpha=0`), incompatible with `alpha=0.110(1)`.  Rosengren is not testable because no complete callable is supplied.
- **Learned:** Universality supplies an independent necessary condition even when low-order coefficient tests already fail.

## E26 — Independent low-temperature series check

- **Question:** Can early LT coefficients be reproduced without finite-lattice Möbius inversion?
- **Method:** A `4x4xc` slab, periodic in the cross-section and `+` on the two transfer faces; first differences in `c` isolate bulk terms via independent exact ghost-bond enumeration.
- **Artifacts:** `experiments/e26_lt_series_independent.py`, `results/series/lt_independent_check.json`, `tests/test_plus_boundary.py`.
- **Outcome:** The independent route agrees with the finite-lattice LT coefficients through `x^15`.
- **Learned:** Agreement across distinct boundary/cluster constructions is more informative than rerunning the same transfer recurrence.

## E27 — Pauli flux centraliser

- **Question:** Can Kitaev-like conserved Pauli fluxes block-diagonalize the natural Ising term generators?
- **Method:** Binary symplectic centraliser proof and brute force on seven graphs.
- **Artifacts:** `experiments/e27_flux_sectors.py`, `results/flux_sectors.json`, `proofs/algebraic_obstruction.md`.
- **Outcome:** For a connected graph, Pauli strings commuting with every individual `X_i` and `Z_iZ_j` are exactly `I` and global spin flip.
- **Learned:** The theorem closes the term-wise Pauli-flux route only.  A later audit corrected an overextension to the full summed-pair commutant; see A01.

## A01 — Hostile replication and adversarial corrections

- **Question:** Which headline claims survive implementations that avoid the primary proof/series code and probe hostile edge cases?
- **Method:** Dense exact matrices, direct torus cycle counts, independent ghost-bond enumeration, benchmark-literal scans, and endpoint recomputation.
- **Artifacts:** `tests/test_audit_replication.py`, `results/audit/audit_replication.json`, `tests/test_plus_boundary.py`.
- **Outcome:** DG formula, tridiagonal residual, direct HT/LT coefficients, and rigorous endpoints replicate.  Three corrections resulted:
  1. `C3` plus an isolated vertex falsifies “DG iff 2-regular”; the proved condition is every degree in `{0,2}`.
  2. Translation on `C3` falsifies “the full commutant has only identity/global parity”; the theorem is only the Pauli centraliser of individual terms.
  3. `plus_boundary=True` undercounted one ghost bond per site when transfer length `c=1`; the implementation was repaired and an independent test now passes.
- **Learned:** Different definitions and edge cases find both implementation bugs and quantifier bugs; counterexamples must be preserved, not silently softened.

## B01 — Bookkeeping acceptance run

- **Question:** What is the actual final test state of the repository snapshot?
- **Method:** `./tests/run_all.sh`, dynamically executing every `tests/test_*.py` with `.venv/bin/python`.
- **Artifacts:** `tests/run_all.sh`, `checkpoints/current_state.md`.
- **Outcome:** At `2026-08-11T19:16:42Z`, **18/18 test scripts passed** in `579.04 s`.
- **Learned:** `test_audit_replication.py` deliberately prints two historical `FAIL` findings while exiting 0: both are successful counterexamples to old prose, and current proof statements carry the corrected scopes.

## 2026-08-13 - Ledger incident and recovery (wave 7)
`notes/hypotheses.csv` was accidentally overwritten (not appended) by the ParametricS agent,
leaving only its own five rows H115-H119 of what had been header + 131 rows. Recovery: the
workspace syncs via Syncthing to the Mac Studio, whose hourly APFS local snapshot
`com.apple.TimeMachine.2026-08-13-062655.local` still contained the full pre-overwrite file
(H001-H164 including all 35 wave-7 rows appended by the seven agents that finished before the
incident; the snapshot predates ParametricS's write, so it contains no H115-H119). The snapshot
was mounted read-only on the Studio, copied back, and the five live H115-H119 rows appended
after it. Validation of the merged ledger: 136 data rows, uniform 11-column CSV, unique
nonempty IDs (checked by `uniq -d`, per-ID grep counts, and a file-level `csv.reader` parse).
Audit backup: `checkpoints/hypotheses_pre_overwrite_snapshot.csv` (removable after final
review). Preventive rule broadcast to all agents: ledger writes are append-only; full-file
rewrites are prohibited.

## 2026-08-13/14 - Waves 7 and 8 (parallel agent rounds; 10 + 18 fronts + audits)
Wave 7 (10 agents): SSOT fix landed first (e29 scoped fields + characteristic_zero_resolution
pointers; canonical artifact regenerated; cross-artifact test validation). Results: full
commutant of (A,B) computed exactly (2x3 = 12, falsifying the dim-2 conjecture; M2+C^8 with a
repeated irreducible and two common eigenlines); support-resolved local-charge theorem with four
genuine 3x3 support-7 charges; interlayer c4 through v^12 (two routes); series second-order/
Mahler/ODE frontiers all negative (138 rows, zero survivors); Kesten ratio inequality honestly
UNRESOLVED (only the limit is in auditable sources) so the lower endpoint stands; KW branches
17/32 empty; Theorem S extended to 3x3; four continuum coupling components (connected interval
honestly unresolved); ladder leading-word induction class exhaustively dead; 2x4 char-0 partial
(sector forms/constituents exact; Killing wall recorded). Ledger incident: ParametricS
overwrote notes/hypotheses.csv; fully restored from the Studio APFS snapshot (see entry above).
Wave 8 (18 agents + 2 audits): headline theorems - complete 2x4 structure
Qz+so(F000)+so(F010)+so(F100)+sl24+sl28+sl4^2 (absolute D21+2B13+A23+A27+2A3, radical=centre,
structural Killing rank 2951) with a clean-room no-producer-import verifier (two blockers fixed
en route: split-form Killing determinant removed as unsupported; per-block central normalization
repaired to one global scale z=(1/8)z_int with exact B-z membership); 2x4 Onsager-quotient no-go
via the exact D21 quotient + Date-Roan; positive-real Kac-Ward cone empty (F=48+8*sum M with a
rational Positivstellensatz identity, all 32 strata); rational KW map 31/32 empty (thin-box
constants), 11111 unresolved with exact negative certificates; ungraded tetrahedron no-go for
arbitrary 8x8 R (complementary minors + Bezout clearance); HT v^24 = 2135670379057/8 exactly
(cut-capped parity frontier, 352 MB vs 242 GB projection, five-prime CRT); interlayer c6 through
v^12 + all-order structure theory (evenness with sharpness witnesses, composition selection
rule, [v^2 w^(2m)]=2 identity, blind cross-agent c6 agreement); commutant family types through
2x5 with 2x6 in [57,71]; pointwise-stabilizer charge criterion (2682-orbit converse, zero
mismatches); infrared barrier theorems (I3/2 is the exact method-class optimum; the K_c
upper endpoint provably method-optimal, the lower endpoint carrying finite obstructions but
not an optimality proof); Lee-Yang 2D calibration passes and 3D non-identifiability
is a certified barrier (sigma=1/2 vs -1/4 both fit exactly); ladder AB/BB family exactly refuted
at L=8 (rank 248/256, eight integer relations); dim_Q(3x2)=263 exact, 3x4 depth-15 >= 1242.
Internal falsification: the proposed degree-four six-mode spectral trace identity was refuted
in-house (three-mode-only necessity; counterexample F(64,1,64,1)=14676480) and converted to a
regression-gated negative artifact before any publication. Two agents (KillingRank24,
SpectralIntervals) died on provider quota after 5 h; e58 remains as a documented unused direct-
Killing route with no artifact; the e76 sector route superseded it. All accepted results were
re-verified by the parent via their standalone tests before ledger append (H165-H259 blocks used
as allocated; 226 rows, unique, 11-column).

## 2026-08-15 - Wave 9 (ten parallel agents, all verified)
Ten fronts from the ranked wave-8 list. Every deliverable standalone-verified by the lead
before its ledger block was appended (H260-H309; ledger now 281 unique rows). Headlines:
3x3 layer algebra completely solved over Q (dim 8034, six sl factors, rad=centre); HT series
to v^26 with external witness and a real square-free-log bug found and fixed; absence intervals
proved unable to span the certified forced-value crossings (real algebraic points), with
three connected crossing-free absence intervals proved; 2x6 commutant closed at 71 = M_7+M_2^2+C^14; first-backtrack SAW
theorem c_(m+n) <= c_m(c_n - a_(n-1)); L=8 ladder defect repaired (D_16 >= 256, patched family)
and the stationary two-letter class proved dead; cabled 16x16 tetrahedron generically rigid
over Q(q); Kac-Ward 11111 got three exact negative certificates (deg<=2 Nullstellensatz
impossible, 16-dim torus reduction, Q(i) mu4 chart empty) but stays open; c6 second route
through v^12 + c4 closed form in 2D data; exact Lee-Yang 5x4x4/5x5x4 zeros with unit-circle
certification, fixed-edge correction class excluded, 5^3 wall quantified. K_c interval
unchanged - wider Simon-Lieb family exactly rejects improvement; missing lemma isolated
(mu^2 <= c_32/c_30 would improve it). No fabricated success; all UNRESOLVED items carry
explicit certificates. Full suite rerun after integration (result in current_state.md).

## 2026-08-15 - Suite-exposed artifact-coupling bug fixed (byte pin -> content pin)
First wave-9 full suite: 74/75; `tests/test_levi_images_2x4.py` failed because it pinned a
byte-level sha256 of `results/algebra_structure/oa_quotient.json`, which `tests/test_oa_quotient.py`
regenerates (fresh `meta` timestamps) on every suite run. The mathematical content was identical;
the pin was fragile by construction. Fix: `e59` producer and its test now pin a canonical
content digest (sha256 of the `data` section, sorted keys); field renamed
`classification_artifact_sha256` -> `classification_artifact_content_sha256`; artifact
regenerated by the producer (11 s). Regeneration-stability proved by running the OA test
(rewrites the artifact) followed by the levi test: PASS. Audit of the same pattern elsewhere:
e62/e86 digests are self-computed or external-immutable - not fragile. Full suite rerun below.

## 2026-08-16 - Wave 10 (ten parallel fronts plus OA synthesis)

Every accepted deliverable was rerun through its standalone verifier by the lead before ledger
append (H310-H361; 333 unique, uniform 11-column rows). Main exact advances: HT free energy to
`v^28` with complete profile CRT certification; the all-`n` SAW theorem
`mu^n<=c_n-a_(n-1)` plus an explicit height-slab `a_35` family, strictly raising the rigorous
`K_c` lower endpoint; exact local uniqueness of all four known spectral forced-value crossings;
complete characteristic-zero rank-drop locus `{ -1,0,1 }` for the canonical cabled tetrahedron;
a Bell(6)/XOR 13-term reduction of `c6` to 2D `G`, connected `U4`, and connected `W6`; the
`3x3` Onsager-quotient corollary; and the finite ladder certificate `D_18>=512`.

Equally important negative closures: the audited finite-torus two-point infrared/GKS LP floor is
positive at every finite even side but `O(1/L)`, so it cannot furnish a uniform magnetization
floor; Kac--Ward branch `11111` has a nonsingular `Q_5` construction component, making rational
Nullstellensatz emptiness impossible at every degree, while that local component fails the
independent order-eight holdout; open `4x4` algebra closure reached only the certified lower bound
1794 and was not promoted to an exact dimension; and the semilinear Lee--Yang evaluator reduced
the `5x5` layer to 2,105,872 representatives but did not produce or claim an isotropic `5^3`
zero. Full-suite acceptance for wave 10: `FINAL: 86 total, 86 passed, 0 failed` (3 h 56 min,
run on a frozen tree while wave 11 worked in isolated worktrees).

## 2026-08-16 - Wave 11 (four isolated worktree fronts, reviewed then re-verified in main)

Method change: each front ran in its own `git worktree` with its own `.venv`, so the wave-10
acceptance suite could keep running on a frozen tree. Every deliverable was then re-verified by
the lead inside the main checkout and reviewed by an independent read-only reviewer before
integration (ledger H362-H365).

Exact advances. (1) **Theorem LC**, the strongest structural result of the program so far in
operator class: on infinite `Z^3`, for every `a,b` with `ab != 0`, the only compact finite-support
operator commuting with `a sum_x X_x + b sum_(<xy>) Z_xZ_y` is a scalar. The proof needs no
computation - a maximal-first-coordinate support site `v` has a private outside neighbour
`w = v+e_1`, whose `Z_w` tensor slot is exactly `b[Q,Z_v]`, after which the `Y_v` slot is exactly
`2ia Q_Z` and the support peels. This supersedes the wave-1 Pauli-string flux theorem on three
axes at once (arbitrary complex observables, a single nondegenerate generator, infinite volume).
(2) **Theorem P**: an ordering-free spectral no-go for the `2x3` layer at `t=1/3`. Counting
collisions among the 2016 unordered eigenvalue-slot pair products replaces inertia windows
entirely; a Gaussian would force gcd degree at least 1351, the monic integral model makes modular
gcd degrees upper bounds, and the exact degree is 385 at two primes. (3) The last high-temperature
Euler-ODE survivor is refuted by its first unobservable coefficient `v^24`, so the HT frontier now
has no unexplained survivor. (4) **Theorem SAW-U**: exact three-family inclusion-exclusion
(`H`, `L`, `C`) gives `a_35 >= 76401062626946721319` and raises the certified `K_c` floor to
`0.2122120589214465859334619330429416874541`.

Two parent audits changed outcomes mid-flight. The dispatch brief for the SAW front asserted that
the height-slab and defect families were disjoint; the lead identified that as unproved (both
families are height-positive), and the corrected certificate proves disjointness from the forced
slab height `h_11 = 7` versus `11-2d`, then goes further by replacing the conservative subfamily
`E` with the full coordinate-monotone family `C` (a gain of exactly `4*2^33`). For the spectral
front, reducing rational coefficients at a prime is not enough - the primitive integral
representative of a monic rational gcd can lose degree - so the certificate was rebuilt on the
integral model `R_int = D R` with `D = 3^15 5^4`, where Gauss's lemma makes the bound hold at
every prime. One review finding (a false "coordinate `v_1+1`" phrasing in generated lemma
metadata, while the proof used the correct "at least `v_1+1`") was fixed at the source, the
artifact regenerated, and a geometric check added to the verifier so the defect class cannot
recur silently.

## 2026-08-16 - Wave 11 acceptance and wave 12 (five isolated fronts)

Wave-11 full suite with the four new tests: `FINAL: 90 total, 90 passed, 0 failed` (2 h 49 min).

Wave 12 kept the isolated-worktree pattern (`/private/tmp/ising3d-w12-*`) and added a standing rule
that every front is reviewed by an independent read-only reviewer *and* re-verified by the lead in
the main checkout before integration (ledger H366-H371).

The wave's central methodological gain was a two-line piece of algebra applied to the wave-11
pair-product invariant. Because the pair-product polynomial is monic in `z` with coefficients in
`Z[t]` after clearing denominators, (i) specialization can only *raise* the derivative-gcd degree,
so a single evaluated coupling bounds the generic value, and (ii) `Res_z(A/G, A'/G)` is a nonzero
element of `Q[t]`, so the couplings where the generic value is not attained form a finite set.
That converted two point certificates into coupling-generic theorems: the `2x3` layer is not a full
six-mode Gaussian spectrum for every `t>0` outside at most `422472960` values, and - the more
interesting one - the physical parity loophole closes for all `t` outside at most `25534083` values,
replacing the wave-6 inertia-based Theorem S' with an ordering-free, coupling-generic statement.
The same invariant also scaled to the `2x4` layer at `t=1/3` (gcd degree `9329` against floor
`26335`), where the inertia route had needed certified symmetry sectors; `3x3` hit a predeclared
90 s wall and carries no claim.

The bound front found a genuinely better construction rather than a bigger search: varying the two
slab rises over `{3,5,7,9,11}` gives 25 slab families that are pairwise disjoint because the heights
after 11 and 22 steps identify the schedule. With two new two-defect gadget languages, exact
inclusion-exclusion gives `a_35 >= 499330428831189067251`, `6.5` times the wave-11 union, and the
certified `K_c` floor rises to `0.2122129322754723621039731646196408804906` (a gain of `8.7e-7`,
about 450 times the wave-11 gain). The lead independently recomputed the block profile with an own
DFS and reproduced the grid total, the grid net, `a_35`, and `M` before integration.

The extensive-charge front answered the question Theorem LC deliberately left open. For
translation-covariant sums of a finite density, vanishing of the derivation is exactly a
lattice-divergence condition, giving a finite exact linear system whose nontrivial quotient is
`rank S_R - rank M_R - 1`. The 1D chain reproduces a growing tower (quotient dimension exactly `2R`
for `R = 0..9`, energy current at `R=1`), while `Z^2` at `R=1,2` and `Z^3` at `R=1` give exactly `1`:
only the Hamiltonian density. That is the first quantitative statement in this program about the
class that actually carries 1D integrability, and it is honestly bounded - `Z^3` at `R>=2` needs
`4^27` density coefficients and was never launched.

Reviews changed three deliverables. The series-grid verifier initially sampled four rows while the
claim quantified over 633 cells; it now replays every cell, every survivor certificate, and all 18
legacy LT rows. The parity certificate asserted that both quotients `A/G` and `B/G` were monic when
`B/G` has leading coefficient `496`; the wording was corrected at the source and is now checked. The
pair-product-scale verifier read the new `2x4` value from the artifact under review instead of
rebuilding it, and the extensive-charge verifier never exercised the rooted-tree current
presentation or its own `R=1` control; both were repaired and re-verified.

## 2026-08-17 - Wave 12 acceptance and wave 13 (eight fronts, first all-size theorem)

Wave-12 full suite: `FINAL: 95 total, 94 passed, 1 failed`. The single failure is
`tests/test_ladder_all_l.py` reporting its own `TimeoutError: NON-DECISIVE: exceeded 3600s` — a
load-induced wall, not a mathematical regression: the machine was simultaneously running eight
research agents, several of the user's unrelated long-lived experiment processes, and my own
verification runs. Rerun in isolation immediately afterwards, the same test PASSES in 1565 s. That
is recorded here rather than papered over: the acceptance number for wave 12 is 94/95 with an
isolated-rerun PASS for the missing case.

Wave 13 ran five computational fronts plus three deliberately high-ambition theory fronts, all in
isolated worktrees, all reviewed and then re-verified by the lead in the main checkout
(ledger H372-H379).

**The headline is Theorem F, the first result in this program with no size restriction.** For every
finite graph with a vertex of degree at least three and every `n`, at inhomogeneous parameters off
an explicit hypersurface, the layer spectrum has at least `3^n - 2^n + 48*3^(n-4)` distinct
unordered-slot pair products and therefore is not a full `n`-mode Gaussian subset-product multiset;
the excess is a constant fraction `16/27` of `3^n`, and Corollary F1 turns it into a
measure-theoretic statement: almost every random-field/random-bond branching layer is spectrally
non-Gaussian, at every size. Two ingredients do the work — an exact tensor identity for the
pair-product count under decoupling, and a single 16-dimensional claw certificate as the seed. The
scope is stated inside the theorem: the physical isotropic curve is not covered, and §9 of the note
*proves* that is an obstruction of the method rather than a gap in the write-up.

A parent hint mattered here. The front first reported that its localization collapsed at uniform
fields; the diagnosis was that `m` decoupled free sites at equal fields have `m` coinciding modes.
Replacing them with an open chain of length `m`, whose modes are pairwise distinct, restored the
exponential excess and produced Theorem U — an isotropic no-go with zero anisotropy at `n = 5..9` —
and then Theorem G, which upgrades all eight isotropic certificates (including `2x3`, `2x4`, `3x3`)
from named points to all-but-finitely-many couplings. Theorem G's `2x3` degree bound came out as
`422472960`, exactly the independently derived wave-12 number.

Adversarial review was decisive on this front, and its findings are worth recording. Theorem U2 and
the claim that multiplicative independence was certified were **withdrawn**: the argument had
substituted modular *lower* bounds into an *upper* bound, which proves nothing. Theorem U's
quantifier was over-general. A "constant fraction of the ceiling" claim used the wrong denominator,
one lemma asserted that `A'/G` is monic when its leading coefficient is `deg A`, and a step in the
1D control placed a nontrivial positive matrix in `SO(2n,R)`, where only the identity is positive
definite. All nine items were corrected, and Theorem F, its corollaries, Lemma 1, Lemma 5, the claw
seed and the `3x3` agreement survived untouched.

Other fronts: the `3x3` layer — the first with a degree-4 vertex — is now certified by a blocked
trace pipeline (gcd degree `59641` against floor `111645`), closing a wall the repository had
recorded as unresolved, with two independent implementations agreeing on the number; the certified
`K_c` floor rose again to `0.2122129498516328253630821339504383555437` via variable-cut slab grids
with 820 explicit pairwise separators; the low-temperature series jumped from `x^32` to `x^52`,
which **refuted 11 of the 18 legacy truncation-unobservable rows** and left the last 7 needing
`x^54`; an exact support-size stratification broke the `Z^3` radius-2 wall for named classes
(quotient exactly 1 for words of at most four non-identity sites); interlayer `c8` landed exact
through `v^12`; and the `2xL` ladder growth conjecture was reduced to a single named saturation
hypothesis whose conclusion, `dim g_L >= 2^(2L-3)`, is quadratically stronger than the conjecture.

Two honest negatives. The all-radius local-charge front found that the general `d >= 2` theorem is
already in the literature (Chiba, arXiv:2412.18903, PRB 111 195130 (2025)); its contribution is an
independent finite-shape reconstruction plus a quantified residual, and a system advisory correctly
caught the first draft of that note stating four finite-shape results as all-radius theorems — the
tags were corrected to `[CONJECTURE]`. And the ladder front's unconditional claim was not retained:
`L = 9` stopped at rank `420/512`, recorded as a non-certificate rather than as evidence either way.

Reviews changed six of the eight deliverables, in every case because a verifier trusted a stored
field or a claim quantified beyond its evidence. That pattern, not any single theorem, is the
reason the numbers in this repository can be trusted.

**Ledger accounting, stated precisely.** `notes/hypotheses.csv` now holds **351 data rows with 351
unique IDs** (verified as `len(rows) == len({r[0] for r in rows})` on `csv.reader` output, not by a
line scan), all uniform 11-column, highest ID `H379`. Row count and highest ID do not coincide
because 28 H-numbers inside earlier per-agent ID blocks were allocated but never used
(`H036, H037, H040, H042, H043, H045..H049, H053..H055, H060, H061, H067..H069, H165..H169,
H190..H194`). There are no duplicate IDs and no embedded newlines: the file has 354 physical lines
= one comment + one header + 351 rows + one trailing blank. An earlier handoff sentence of mine
said "351 rows, 353 IDs", which is arithmetically impossible with one ID per row; this is the
corrected statement.

**Wave-13 acceptance suite verdict.** `FINAL: 103 total, 102 passed, 1 failed` (14 h 25 m). The
single failure is `tests/test_ladder_all_l.py` reporting its own `TimeoutError: NON-DECISIVE:
exceeded 3600s` — the wall-clock gate measuring the machine's load average (85–118 at the time),
not the computation. Rerun in isolation the same test PASSES in 1565 s. This is the same resource
artefact that already appeared in the wave-12 run; the deterministic fix follows below.

**Deterministic compute gating, applied.** Every per-test compute budget is now CPU time
(`time.process_time()`), never wall clock: `tests/test_ladder_all_l.py` replaces its
`signal.alarm(3600)` wall gate with an audited `process_time()` deadline (the SIGALRM remains only
as a 4x safety net), `tests/test_pair_product_3x3.py` and `tests/test_extensive_3d_r2.py` replace
their `perf_counter()`/`signal.alarm` budget checks the same way, and `tests/test_interlayer_c8.py`
replaces its `monotonic()` guard. All four converted tests then PASS in the main checkout on a
loaded machine: spectral 541 s CPU, interlayer 228 s CPU, extensive 981 s wall (its gate is
informational), and the ladder test by construction can no longer fire on elapsed load. RSS limits
stay wall-independent (peak memory is not a function of other processes' presence). The
"recheck in isolation" convention remains only as reading guidance for older runs.

**Wave-13 verification evidence, per test, all run by the lead in the main checkout with the
repository interpreter.** Recorded here because the 103-script acceptance suite (`fullsuite13`) runs
for many hours under the machine's unrelated background load, and several of these tests carry their
own wall budgets:

| test | result | wall |
|---|---|---|
| `tests/test_pair_product_3x3.py` | PASS | 596 s |
| `tests/test_saw_union3.py` | OK (25 checks) | 131 s |
| `tests/test_lt_x34.py` | OK | 72 s |
| `tests/test_interlayer_c8.py` | PASS | 257 s |
| `tests/test_extensive_3d_r2.py` | PASS | 981 s |
| `tests/test_allr_charges.py` | OK (19,584 assertions) | 241 s |
| `tests/test_allsize_gaussian.py` | 53 passed, 0 failed | 165 s |
| `tests/test_ladder_alll_proof.py` | PASS (178 checks) | 2,808 s |

Independent arithmetic re-derivation by the lead, from the stored artifacts only:
`a_35 = 507841195880595165555`, `M = c_36 - a_35 = 2940863015138821131395115`,
`K_c >= 0.212212949851632825363082133950438355543751098`; Theorem F's excess ratio
`48*3^(n-4)/3^n = 0.592593 = 16/27` at `n = 4, 6, 9, 12`; and the `3x3` thresholds
`C(512,2) = 130816`, `3^9-2^9 = 19171`, floor `111645`, `130816 - 59641 = 71175`.
