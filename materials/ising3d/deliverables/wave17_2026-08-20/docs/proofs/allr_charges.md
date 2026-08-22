# Structure of translation-covariant extensive charges: finite-shape certificates

Artifacts: `experiments/e127_allr_charges.py`, `tests/test_allr_charges.py`,
`results/integrability/allr_charges.json`.

Reproduce:

```
PYTHONPATH=src .venv/bin/python experiments/e127_allr_charges.py
PYTHONPATH=src .venv/bin/python tests/test_allr_charges.py
```

**Status convention, and a warning about it.** `[COMPUTATION]` is an exact finite
calculation recorded in the artifact and independently rebuilt by the standalone
verifier (19,584 assertions). `[EXTERNAL]` is an audited outside result.
`[CONJECTURE]` is a general statement that the finite data supports but that is
**not proved here**. `[UNRESOLVED]` marks what is deliberately not decided.
Everything below that quantifies over *all* radii or *all* shapes is either
`[EXTERNAL]` or `[CONJECTURE]`: this note contains no all-radius proof of its
own. An earlier draft of this file overstated four items as theorems; the
statements are unchanged, the tags are corrected.

---

## 0. Attribution, stated first

**[EXTERNAL]** The general absence of local conserved quantities for this model
class in `d >= 2` is **not new here**. It is due to Y. Chiba,
*Absence of local conserved quantities in the transverse-field Ising and related
models in two or more dimensions* (arXiv:2412.18903; Phys. Rev. B **111**,
195130 (2025)). Chiba works on a periodic lattice of linear size `L` with the
support-radius hypothesis `k_mu* <= L/2` and position-dependent coefficients; his
Propositions 1 and 2 carry the argument. That reference, not this note, is the
proof of the unrestricted `d >= 2` statement.

**What this file contains.** An independent finite-shape reconstruction of the
mechanism in the `Z^d` translation-orbit framework used elsewhere in this
repository, machine-checked shape by shape; exact single-shape kernel
computations in `d >= 2` and their `d = 1` counterexamples; and an exact
quantification of the residual kernel that the leading-symbol argument alone does
not remove. No priority and no general proof is claimed.

---

## 1. Setting

For `H = field * sum_x X_x + bond * sum_(<xy>) Z_x Z_y` on `Z^d` and a finite
density `q`, the formal translation sum `Q(q) = sum_x tau_x(q)` is conserved iff
`pi(1/2)[H,q] = 0`, where `pi` is the translation-orbit sum; this divergence
criterion is the one proved in `proofs/extensive_charges.md`. Write
`Lbar = pi . (1/2)[H, .]`. For a word `X^a Z^b` let its *bounding box* have shape
`sigma = (sigma_1, ..., sigma_d)`, let `V_sigma` be the span of words of shape
exactly `sigma`, and for direction `i` call the two extreme hyperplanes of the box
the *extreme `i`-faces*.

## 2. Extent grading and the two-entry structure

**[COMPUTATION — extent grading]** On every tested shape, the shape-`(sigma+e_i)`
component of `Lbar` restricted to `V_sigma` coincides with the direction-`i`
growth part of the commutator, checked term by term
(`leading_symbol_identification_lemma_B`, `all_images_in_target_shape`).
Structural reason, valid for a single word: the field term
`(1/2)[X_s, X^a Z^b] = [b_s=1] X^(a+e_s) Z^b` cannot enlarge the box because
`b_s = 1` already puts `s` in the support, and the bond term
`(1/2)[Z_u Z_v, X^a Z^b] = -[a_u xor a_v = 1] X^a Z^(b+e_u+e_v)` enlarges it only
when exactly one endpoint lies outside, and then by one layer in one direction.
This is a per-word observation; the orbit-level identification is what the
artifact checks, and it is checked only on the recorded shapes.

**[COMPUTATION — two nonzero entries]** Every orbit row of that leading symbol has
at most **two** nonzero entries (`max_nonzero_entries_per_row = 2` on all 15
recorded shapes). Its kernel is therefore a grounding/ratio-propagation problem
on a graph whose components are grounded (forced zero) or free, and an
independent Gaussian-elimination solver reproduces the graph solver exactly
(`graph_solver_matches_gaussian_elimination`).

## 3. The two forcing statements

**[CONJECTURE T1]** In a `sigma_i`-maximal layer, any word whose `X`-support meets
**both** extreme `i`-faces has coefficient zero.
**[COMPUTATION]** Verified on every recorded shape, with `T1_words` counted and
`T1_all_forced_zero` true (`T1_holds_on_every_shape`).

**[CONJECTURE T2 — the `d >= 2` clause]** A word whose `X`-support meets exactly
one extreme `i`-face, and whose **other** extreme `i`-face carries at least two
support sites, has coefficient zero.
**[COMPUTATION]** Verified on every recorded shape (`T2_holds_on_every_shape`).

**[COMPUTATION — T2 is where `d >= 2` enters]** On a chain every face is a single
site, so T2's hypothesis is vacuous: `T2_words = 0` for every chain shape
`(2) .. (10)` and `T2_words > 0` for every `d >= 2` shape
(`d1_control_T2_hypothesis_vacuous`, `dge2_T2_hypothesis_has_content`). This is
the discriminating check — a mechanism that also applied in `d = 1` would be
wrong, because the chain genuinely has a charge tower.

**[CONJECTURE M]** After T1 and T2, the surviving words of the maximal layer are
exactly the interior-`X` words together with the thin-face words whose single
boundary site carries `Z` with an `X` on its inward neighbour.
**[COMPUTATION]** Verified pointwise, row by row, on every recorded shape
(`local_rows_verify_theorem_M_pointwise`), with the `d=1` / `d >= 2` controls
separated (`local_d1_control_no_thick_faces`,
`local_d1_control_clause_b_words_present`, `local_dge2_thick_faces_dominate`).

## 4. Single-shape densities: exact finite kernels

**[CONJECTURE S]** A density all of whose words share one bounding-box shape
carries no extensive charge in `d >= 2`.

**[COMPUTATION]** The exact kernel of the combined leading-symbol system,
restricted to that class, is `0` at shapes `(2,2)`, `(3,2)`, `(2,3)`, `(4,2)`,
`(2,4)`, `(3,3)`, `(2,2,2)`, at two coupling pairs each
(`dge2_no_single_shape_charge_density`). These seven shapes are the certificate;
they are not an induction, and no general argument is supplied here.

**[COMPUTATION — the mandatory `d = 1` failure]** The same code on chains gives
kernel dimension exactly **1** at every shape `(2) .. (7)`, with the explicit
conserved witness family

```
    j_n = Y_0 X_1 ... X_(n-2) Z_(n-1)  -  Z_0 X_1 ... X_(n-2) Y_(n-1),
```

verified conserved for `n = 2 .. 8` (`chain_single_shape_charge_family_conserved`,
`d1_single_shape_charge_kernel_is_one_dimensional`). So Conjecture S is false in
`d = 1`, exactly as it must be.

## 5. Global quotients, cross-checked

**[COMPUTATION]** Independently recomputed nontrivial quotient dimensions agree
with the wave-12 census (`global_quotients_are_2R_in_d1_and_1_for_d_ge_2`,
`cross_reference_e117_quotients`): `d = 1` gives `0, 2, 4, 6, 8, 10` at
`R = 0..5` — the chain tower — while `d = 2` gives `0, 1` at `R = 0, 1` and
`d = 3` gives `0` at `R = 0`.

## 6. The residual obstruction, quantified

**[COMPUTATION]** The all-direction leading symbol alone does **not** close the
general `d >= 2` case (`leading_symbol_alone_does_not_close`). Its active kernel
remains nonzero:

| shape | residual active kernel dimension |
|---|---:|
| `(2,2)`, `(2,2,2)` | 2 |
| `(3,2)`, `(2,3)` | 6 |
| `(4,2)`, `(2,4)` | 28 |
| `(3,3)` | 184 |
| `(5,2)`, `(2,5)` | 120 |

These residual words are the analogues of Chiba's type ii-b / iii-b terms.
Removing them requires uniform control of the **shape-preserving subleading**
level, which is exactly what this note does not have; that gap is why T1, T2, M,
and S are tagged `[CONJECTURE]` rather than proved here.

## 7. Scope

**[UNRESOLVED]** Nothing here classifies quasilocal charges, non-translation-
covariant charges, densities outside the tested shapes, or non-linear-in-`H`
charges, and nothing here proves an all-radius statement. No non-integrability
theorem for the 3D Ising model follows, and no thermodynamic-limit or spectral
claim is made. For the unrestricted `d >= 2` theorem, cite Chiba (sec. 0).
