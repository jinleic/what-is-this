# Adversarial audit that changed the repository: five claims corrected, one code bug fixed

## Referee verdict

**The repository does not solve the three-dimensional Ising model.** Its own final report now says
so. Several finite-volume theorems and exact computations are valuable, but none produces a
thermodynamic-limit free energy, a closed critical manifold, or exact 3D critical exponents.

This audit did not merely endorse the repository. It found **six material defects**:

1. **A false corollary of Theorem DG.** The defect does not vanish only on 2-regular graphs: an
   isolated vertex also contributes zero. The exact counterexample `C3` plus one isolated vertex
   has a zero dense-matrix defect. The statement and proof were corrected to `deg(v) in {0,2}`.
2. **A published benchmark entered a fit-selection rule.** `experiments/e04_lee_yang.py` originally
   selected “high-temperature” inputs using `K < 0.221654626` and then fit edge data. That is an
   input to the analysis, despite not being an adjustable parameter. It was replaced by the
   repository's independently certified lower bound on `K_c`; the Lee–Yang checks still pass.
3. **A real boundary-condition code bug.** `box_broken_bond_poly(..., plus_boundary=True)` omitted
   one of the two outward transfer-direction ghost faces when the transfer length was one. The
   `(2,2,1)` polynomial failed an independent spin enumeration, first at `x^4`. The source was
   repaired to apply both ghost faces to the unique layer. The current implementation agrees with
   independent enumeration on all 22 cases in `tests/test_plus_boundary.py`; the FLM coefficients
   are unchanged because shape canonicalisation had isolated this case.
4. **The Pauli-centraliser theorem was given an invalid full-commutant consequence.** On `C3`, the
   one-site translation `U` commutes exactly with both summed generators, satisfies `U^3=I`, and
   `{I, prod X, U}` has exact Gram determinant `480`. Thus the full commutant is not two-element and
   the Hilbert space does not have only two symmetry sectors. The text was narrowed: Theorem 11
   excludes same-Hilbert-space **Pauli-string** flux families of the Kitaev kind, not space-group
   sectors, non-Pauli fluxes, or ancilla-extended mechanisms.
5. **A scalar Kac–Ward no-go was over-quantified.** The proof established a contradiction only
   for the strict two-parameter cubic-covariant scalar rule, while one sentence said it covered
   “any scalar turn rule.” The full 30-weight scalar family was explicitly unresolved. The text
   now carries the correct strict-covariance qualifier and records the unresolved family.
6. **A final synthesis falsely merged independent obstructions.** It said all four no-go results
   pointed to the quartic `YZZZ` defect, but the Pauli-flux theorem is degree-independent and also
   holds in 1D where that quartic term is zero. The final report now says only DG, TD and the
   claw/Gaussian obstruction share the quartic object; Theorem 11 closes a distinct route.

The most important positive result survived: a completely separate dense-matrix implementation
reproduces the closed Dolan–Grady defect exactly on `K1,3` and the open `2x3` grid. The dense
tridiagonal test on the `3x3` torus also gives a large, exactly nonzero least-squares residual.

Machine-readable evidence is in `results/audit/audit_replication.json`, produced by
`tests/test_audit_replication.py`.

## Tags and evidentiary standard

* **T — theorem:** a mathematical derivation is supplied locally.
* **M — machine-verified:** exact integer/rational computation unless explicitly labelled
  floating point.
* **N — numerical:** a floating-point observation or finite-series inference.
* **E — external:** imported from cited literature; this audit did not turn it into a local proof.
* **C — conjecture:** extrapolation or unproved structural expectation.

“Survives” below means that the particular attack stated in this report did not break the claim. It
does not promote a finite computation to a theorem. Context-designated baseline checks (for example
the 17 enumeration/transfer comparisons) are cited rather than wastefully rerun.

## Complete claim inventory

The inventory groups repeated JSON renditions of the same mathematical statement, but names every
claim-bearing file family.

| claim | tag | source(s) | hostile disposition |
|---|---:|---|---|
| Model normalisation `Z=e^{Kn_b} sum_q c_q x^q = 2^N cosh(K)^{n_b} P(v)` and the length-1/length-2 bond conventions | T/M | `problem_specification.md`, `notes/CODE_API.md`, lattice and transfer modules | **Survives current tests.** Independent bond multisets and broken-bond polynomials agree on degenerate tori. One separate plus-face bug was found and repaired. |
| Exact spin enumeration and exact even-subgraph polynomial | M | exact-enumeration module; `tests/test_tm_vs_enumeration.py` | Context-established exact baseline; 17 cross-implementation cases. Finite-volume only. |
| Exact layer transfer matrix | M | transfer module; `tests/test_tm_vs_enumeration.py` | Context-established baseline. The public plus-boundary path had the repaired `c=1` defect. |
| Dimensional-reduction, zero/infinite-coupling, degeneracy, spin-reversal and cycle-space gates | M/T | `tests/test_gates.py`, final report | Finite identities; no thermodynamic 3D solution follows from them. |
| Kaufman four-sector finite-torus formula | N/M/E | `notes/onsager_derivation.md`, `results/onsager/onsager_2d.json` | Reproduces exact integers to about `3.6e-60`; original paper formula was unobtained and the sector convention was inferred against exact data, so “independent reproduction” is weaker than source-level derivation. |
| Onsager double-integral free energy and `sinh(2Kx)sinh(2Ky)=1` | T/N/E | same files | The thermodynamic integral, Catalan critical value and `K_c=asinh(1)/2` controls survive. The literal `0.4406867...` occurs in finite-size control inputs, not the derivation of `K_c`. |
| Free-fermion/Majorana form of the 2D transfer operator | N/E | Onsager files | Off-bilinear residual near `1.2e-15`; numerical Clifford diagnostic, not an exact symbolic identity. |
| Lee–Yang polynomial palindrome, root count and circle location on 13 finite lattices | M/N | `notes/lee_yang.md`, `results/lee_yang/lee_yang_analysis.json` | Exact coefficient symmetries survive; root radii are high-precision numerical results. `tests/test_lee_yang.py` was rerun after benchmark decontamination and passed. |
| Lee–Yang circle theorem for all ferromagnetic volumes | E | Lee–Yang citations | The finite checks do not prove the external theorem; the theorem supplies the general statement. |
| First-zero scaling and Yang–Lee edge fits | N | Lee–Yang artifact | Exploratory finite-size fits only. The original benchmark-dependent sample selection was invalid as an “input-free” analysis and was repaired. No exact exponent is established. |
| Exact finite-volume 3D Ising–`Z2` gauge duality, including twists and orbit factors | T/M/E | `proofs/duality_3d.md`, `results/duality/duality.json` | 27 exact checks and the chain-complex derivation survive for the tested boxes. The 3D cases have all side lengths at least two; periodic length one is not covered. |
| “3D Ising is self-dual” | falsified T/M | duality proof/artifact | Correctly rejected: the dual model is a one-form gauge theory and the torus requires sector sums. |
| Thermodynamic Ising/gauge free-energy relation | T/E | end of `proofs/duality_3d.md` | Conditional on existence and boundary/twist-sector equivalence in the limit. The exact finite identity alone is not a 3D critical-point equation. |
| Finite-lattice Möbius inversion for `log P` / `log Xi` | T/M | `notes/flm_derivation.md`, series module | Algebraic identity survives. Exactness order depends entirely on the two geometric bounds audited below. |
| HT order bound `2[(a-1)+(b-1)+(c-1)]` | T | FLM note/artifacts | **Survives an independent proof**, including Mayer clusters, given below. |
| LT order bound `4(a+b+c)-6` | T | FLM note/artifacts | **Survives an independent projection proof**, including Mayer clusters, given below. |
| Simple-cubic HT free energy through `v^20` | M | `results/series/sc_ht_free_energy.json`, `extended_sc_ht_free_energy.json` | Exact Fractions/CRT. Independently spot-checked at `v^2,v^4,v^6`; repository also has a separate GF(2) route through `v^10`. Orders above that remain dependent on the FLM implementation and CRT bookkeeping. |
| Simple-cubic LT free energy through `x^28` | M | `sc_lt_free_energy.json`, `extended_sc_lt_free_energy.json` | Exact Fractions. Independently spot-checked at `x^6,x^10,x^12`; repository slab calculation checks through `x^15`. Higher terms have no second implementation. |
| Series-only `K_c=0.2229 +/- 0.0050` | N | `notes/series_analysis.md`, `extended_series_analysis.json` | Benchmark is used only after the current approximant ensemble is formed. Result is deliberately imprecise and method-dependent, not exact. |
| Series exponent `alpha=0.35 +/- 0.34` | N | same files | Dlog-Padé and differential approximants disagree strongly. The wide error bar honestly exposes that the series is not informative about `alpha`. |
| No low-order Euler-form D-finite ODE | M/C-scope | same artifact | Exact finite search with held-out coefficients and an exponential control. It excludes only the listed `(r,d)` budget and is not non-D-finiteness. |
| SAW lower bound on `K_c` | T/M/E | `proofs/kc_bounds.md`, `results/bounds/kc_bounds.json` | Inequality direction and outward endpoint survive. Uses external correlation domination plus exact `c_14=4468911678`. |
| Square-layer upper bound | T/E | same | Direction survives GKS monotonicity: adding positive interlayer bonds orders earlier, hence `K_c(3D) <= K_c(2D)`. |
| Infrared upper bound `K_c <= W_sc/6` | T/N/E | same | Normalisation and outward endpoint survive independent 110-digit evaluation. The decisive Gaussian-domination theorem and Watson gamma identity remain external inputs. |
| Simon–Lieb finite-box lower bounds | T/M/N/E | `proofs/simon_lieb_bound.md`, `simon_lieb_bounds.json` | Criterion hypotheses and safe direction are stated. Small rational cases are exact; long-prism values use a conservative positive-recurrence forward-error model. The result is weaker than the SAW bound and is not used in the final endpoint. |
| Theorem DG closed formula | T/M | `proofs/dolan_grady_defect.md`, `dg_theorem_check.json`, `dg_local_structure.json` | **Survives independent exact dense matrices** on the two required graphs. |
| Original “DG defect vanishes iff 2-regular” | falsified | historical theorem text | Broken by `C3` plus an isolate; now corrected to degrees in `{0,2}`. |
| Quartic defect count `sum_v C(deg(v),3)` | T/M | DG proof/artifacts | Survives the dense cases and follows from distinct Pauli support. It is an obstruction to the Dolan–Grady mechanism, not a proof of non-solvability. |
| No tridiagonal/Askey–Wilson relation on each tested layer | M/T-finite | `proofs/tridiagonal_nogo.md`, three tridiagonal artifacts | Exact short-string certificates cover the named layers through `5x5`. Dense `3x3` least squares independently survives. |
| No TD relation for every possible 2D layer size | C/not claimed after scope correction | TD proof/final report | **Not proved.** The repository now explicitly limits the statement to tested sizes. |
| Claw obstruction and subdivision theorem | T/M/E | `proofs/algebraic_obstruction.md`, `dla_table.json` | Self-contained necessity result for term-wise Majorana bilinearity survives. It concerns the natural term representation, not every factorisation of the transfer operator. |
| DLA dimension bound for an ancilla Majorana representation | T | same proof | Algebra-dimension implication is valid under its stated full conjugation hypothesis. |
| Exact DLA dimensions for finite layers | M | `dla_table.json` | Finite exact symplectic closure. The `3x3` value forces 128 modes under Theorem 7. |
| Exponential DLA / two-generator algebra growth in lattice families | C | `onsager_algebra_dims.json`, `minimal_obstruction.json`, algebra notes | Data are striking but finite: `2x2,2x3,2x4` and a few perturbed paths cannot prove asymptotics. The final report now labels the growth law conjectural/open. |
| Finite one-dimensional algebra formulas `n^2` and `3n-1` | M/C | algebra artifacts | Verified to finite `n`; not a theorem without a proof. |
| Pauli-group centraliser `{I,prod X}` | T/M | `proofs/algebraic_obstruction.md`, `flux_sectors.json` | Elementary symplectic proof survives for connected graphs and Pauli strings. |
| Original “only two invariant sectors/full commutant” consequence | falsified | historical proof/final report | Dense `C3` translation counterexample, Gram determinant `480`; scope now corrected to Pauli-string Kitaev fluxes. |
| Original “all four no-gos point to the quartic defect” synthesis | falsified scope claim | historical final report | The Pauli-flux theorem holds even in 1D where the quartic defect vanishes; corrected so only DG, TD and the claw/Gaussian obstruction are grouped. |
| No local conserved charges at ranges 2–4 on the tested 2D layer | M | `notes/integrability.md`, `conserved_charges.json` | Exact finite-field nullities for that ansatz. It does not exclude longer-range, nonlocal or differently structured charges. |
| Transfer-matrix commutant multiplicities | N | same artifact | High-precision spectral clustering on finite matrices, not an exact commutant classification. |
| Bond-dimension-two tetrahedron-equation failure | M/C-scope | `tetrahedron.json`, integrability note | Exact polynomial specialization rules out the tested tensor/gauge family. Higher bond dimensions, alternating tensors and other vertex formulations remain open. |
| Planar Kac–Ward identity `det(I-v Lambda)=P(v)^2` on free square graphs | E theorem + M control | `proofs/kac_ward_3d_obstruction.md`, `results/kac_ward/kac_ward.json` | Exact modular reconstruction gives zero discrepancy on four boxes; the all-planar-graph statement remains the external Kac–Ward/Whitney theorem. |
| No strict cubic-covariant scalar turn rule in 3D | M | same proof/artifact | Exact rational Gröbner basis becomes `{1}` at order `v^8`; the fitted finite-box candidate independently fails at the same order. |
| No full 30-weight scalar Kac–Ward rule | C/open | same | **Not proved.** The generic gauge-chart system is consistent through order six and timed out at order eight; zero branches, projective weights, spin-structure sums and bond dimension at least two are untested. |
| Original “no scalar turn rule” wording | falsified scope claim | historical Kac–Ward proof text | Corrected from “any scalar turn rule” to “any strict cubic-covariant scalar turn rule.” |
| Zhang 2007 critical point and magnetisation signs | E target + M/N falsification | `notes/falsification.md`, `falsification.json` | Benchmark discrepancy plus an independently exact low-order series refute the claims as transcribed. |
| Degang Zhang 2021 anisotropic critical line/free energy | E target + T/M falsification | same | Axis-permutation failure is decisive; exact HT mismatch at `K^4` is independent support. |
| Rosengren critical-point-only conjecture | E/N | claimed-solution registry/falsification note | Differs from the benchmark, but it never supplied a free energy. This is a numerical falsification target, not an exact derivation. |
| Bootstrap critical dimensions/exponents and amplitude ratios | E | `results/cft/cft_constraints.json`, literature matrix | Imported literature data. Algebraic scaling identities among quantities derived from the same two dimensions are consistency checks, not independent measurements. |
| Claimed-solution formulas fail the 3D `alpha` constraint | N/E | CFT artifact | The extractor passes 2D/synthetic controls and rejects the two available formulas; nevertheless this is a numerical singularity classification of already-falsified candidates. |
| Current series determine universal amplitude ratios | falsified/absent | CFT artifact | Correctly reports none: the zero-field series lack controlled two-sided continuation and field derivatives. |
| An exact 3D free energy or exact critical exponents | C/open | `problem_specification.md`, final report | **Not obtained.** |

## Attacks actually carried out

### 1. Dense Theorem DG replication: survives

I constructed `A` and `B` directly as dense Kronecker-product matrices. No repository Pauli-string,
symplectic or commutator code is imported. Because `iY` is a real integer matrix, every entry and
every matrix product is exact `int64` arithmetic at these sizes.

For both required graphs I computed

```
R = [B,[B,[B,A]]] - 16[B,A]
R_closed = 24 i sum_v Y_v ((deg(v)-2) sum_{u~v} Z_u
                            + 2 sum_{triples in N(v)} Z_u Z_u' Z_u'')
```

and independently checked `[A,[A,[A,B]]]-16[A,B]`.

| graph | matrix size | `max abs(R)` | `max abs(R-R_closed)` | first DG residual |
|---|---:|---:|---:|---:|
| star `K1,3` | `16x16` | 120 | **0** | 0 |
| open `2x3` grid | `64x64` | 120 | **0** | 0 |

This is stronger than a `1e-12` numerical agreement: it is exact zero by a separate dense path.

The same code found the failed corollary: `C3` plus an isolated fourth vertex is not 2-regular but
has `max abs(R)=0`. The corrected theorem's degree condition `{0,2}` is consistent with it.

### 2. Dense Theorem TD replication: survives for the `3x3` torus

For the independently constructed 18-edge `3x3` periodic square graph, I built `512x512` dense
integer matrices and the four TD2 operators

```
T0 = [B, B^2 A + A B^2]
T1 = [B, B A B]
T2 = [B, B A + A B]
T3 = [B, A].
```

The three design columns have rank 3. Solving the normal equations over `Fraction` gives the unique
least-squares minimiser

```
(beta, gamma, rho) = (23/13, 18/13, 508/13).
```

The exact minimum squared Frobenius residual is

```
488374272/13,
```

so the minimum norm is `6129.212974951`. The target norm is `25051.401876941`, hence the relative
residual is `0.244665468426`. The largest residual entry is exactly `2880/13`, with 2880 nonzero
entries after clearing denominators. This is nowhere near a roundoff-level zero. It independently
confirms “no TD2 parameters” for this graph, but says nothing by itself about untested layer sizes.

### 3. Periodic length-one/length-two and prefactor attack: current code survives; one bug was found

`tests/test_audit_replication.py` independently constructs the specified bond multiset and enumerates
all spins for periodic shapes `(1,1)`, `(1,2)`, `(2,2)`, `(1,2,3)` and `(2,2,3)`. It checks both the
lattice bond lists and the transfer polynomial. Length one contributes no self-loop; length two
contributes two parallel bonds. Every current case agrees exactly.

The distinct plus-boundary test initially failed on `(2,2,1)`: the unique transfer layer received
only one outward `z` ghost face. This was a genuine source bug, not a convention preference. After
repair, my oracle agrees coefficient-for-coefficient, and the dedicated independent test passes 14
3D plus 8 2D boxes. The HT/LT series were regenerated and did not change.

The prefactors in the exact enumeration, transfer, FLM and duality notes are consistent with the
problem specification. No thermodynamic conclusion was accepted merely because a finite prefactor
matched.

### 4. Circularity and literal-constant scan: one leak found and repaired

Executable Python under `src/`, `experiments/` and `tests/` was scanned for `0.221654626`,
`0.4406867` and `0.629971`.

* `e04_lee_yang.py` **did** use `REFERENCE_K` to choose inputs to its edge fit. That violated the
  blanket “evaluation only” description. The current code selects only points below
  `0.2074277114992...`, the locally certified lower bound on `K_c`, so every selected point is
  rigorously disordered without importing the benchmark.
* `e19_series_analysis.py` computes its retained approximants, central estimate and stability
  envelope before `benchmark_difference`; the benchmark is a subsequent comparison.
* `e05`, `e06`, `e12` and their tests use the 3D value in bound/falsification comparisons.
* `e02` and `test_onsager.py` include the known 2D critical value as one control evaluation, while
  the critical coupling used in the derivation is independently `asinh(1)/2`.
* No executable occurrence of `0.629971` was found; it enters the CFT artifact as external data
  derived from quoted scaling dimensions.

A literal scan cannot detect a human choosing an analysis after seeing the benchmark. What it can
establish is that the current code has no direct data flow from the 3D benchmark into either series
fit or Lee–Yang sample selection.

### 5. Independent proof of the HT FLM order bound: survives

Let a connected Mayer cluster contributing to `log P` span side lengths `a,b,c`; its coordinate
spans are `a-1,b-1,c-1`. Regard every even-subgraph polymer with its edge multiplicity. Because the
polymer incompatibility graph is connected, their union is a connected Eulerian multigraph (shared
vertices join the Eulerian pieces). It therefore has a closed Euler tour.

Project that closed tour onto coordinate `i`. A closed integer walk that reaches both the minimum
and maximum coordinate must make at least `span_i` positive and `span_i` negative steps. Thus the
number `m_i` of edges in direction `i`, counted with polymer multiplicity, obeys

```
m_i >= 2 span_i.
```

The activity order is the total edge multiplicity, so

```
order = m_x+m_y+m_z
      >= 2[(a-1)+(b-1)+(c-1)].
```

This proof covers the Mayer-cluster case, not merely one connected even graph. Hence boxes whose
minimum possible order exceeds the truncation cannot affect the retained HT coefficients.

### 6. Independent proof of the LT FLM order bound: survives

Let `U` be the connected union of the droplets in a contributing Mayer cluster. The sum of their
surface activities is at least the exterior edge boundary `|partial U|`: every boundary edge of the
union appears in at least one constituent surface.

For the `x` direction, every nonempty `x`-line has at least two boundary edges, so

```
|partial_x U| >= 2 |projection_yz U|.
```

The `yz` projection of a connected set is connected and spans lengths `b,c`; any connected subset
of the square lattice with those spans contains at least `b+c-1` projected sites. Cyclically,

```
|projection_yz U| >= b+c-1,
|projection_xz U| >= a+c-1,
|projection_xy U| >= a+b-1.
```

Summing the three directional boundary estimates gives

```
order >= |partial U|
      >= 2[(b+c-1)+(a+c-1)+(a+b-1)]
      = 4(a+b+c)-6.
```

This also handles overlapping/incompatible formal polymers because surface multiplicity can only
increase relative to the boundary of the union. The claimed LT truncation rule is therefore in the
safe direction.

### 7. Direct series spot checks without FLM: survive

For HT orders through six I used a `7x7x7` torus. A winding cycle of length at most six would need a
period no larger than six, so `L=7` excludes wrap-around contamination. Direct simple-cycle counts
are

```
# 4-cycles = 1029 = 3 per site,
# 6-cycles = 7546 = 22 per site.
```

Together with the exactly expanded `3 log(cosh K)` prefactor, this gives coefficients of
`phi-log 2`

```
[v^2, v^4, v^6] = [3/2, 15/4, 45/2],
```

and interaction coefficients `[v^4,v^6]=[3,22]`, exactly as claimed.

For LT through order 12 I used a `4x4x4` torus. `L=3` would be unsafe because a wrapping three-spin
cycle has surface 12; at `L=4` the shortest wrapping droplet loop has surface 16. Direct surface
counting gives one monomer per site at surface 6, three nearest-neighbour dimers per site at surface
10, and `N+2E=7N` incompatible ordered monomer pairs. Therefore

```
[x^6, x^10, x^12] = [1, 3, -7/2].
```

This code imports neither the transfer matrix nor the series module.

### 8. Rigorous-bound direction and endpoint attack: survives with external dependencies

The SAW direction is correct. Submultiplicativity gives
`mu <= c_n^(1/n)`, hence `c_n^(-1/n) <= 1/mu`. High-temperature correlation domination gives
`1/mu <= v_c`, so

```
atanh(c_14^(-1/14)) <= K_c.
```

Adding ferromagnetic interlayer bonds can only increase correlations, so the simple-cubic model
orders no later in `K` than decoupled square layers: `K_c(3D) <= K_c(2D)`.

For the infrared route, the stated normalisation yields spontaneous order whenever
`K > W_sc/6`, hence `K_c <= W_sc/6`. I independently evaluated the two formulas at 110 digits:

```
SAW point = 0.2074277114992039908436804465100887760027696...
IR point  = 0.2527310098586630030260020266135701299925645...
```

The recorded decimal endpoints are outward, and

```
0.2074277114992039908436804465100887760027
  < 0.221654626
  < 0.2527310098586630030260020266135701299926.
```

No unbounded numerical quadrature is used in the final interval. However, this audit did not
re-prove reflection positivity, the infrared bound, the SAW correlation inequality, or the
Glasser–Zucker Watson-integral identity; those remain external theorem dependencies. The long-box
Simon–Lieb result has an explicit positive-recurrence forward-error allowance, but a future
certificate should replace the operation-count model with exact rationals or outward interval
arithmetic for the final reported prism.

### 9. Finite-volume versus thermodynamic-limit attack: limited claims only

* DG, TD, DLA, conserved-charge, tetrahedron and flux results are finite-layer statements. A
  nonzero defect density is evidence against a mechanism, not a theorem that no exact bulk free
  energy exists.
* The duality identity is exact at finite volume only after the full twist-sector transform and
  gauge multiplicities are included. Passing to a bulk free-energy relation assumes existence and
  boundary/twist independence; it supplies no self-dual scalar equation in 3D.
* Kaufman/Onsager is the one track where a thermodynamic-limit integral is actually taken and
  checked against the independent critical Catalan form.
* FLM coefficients are thermodynamic coefficients because the geometric order bounds prove finite
  stabilisation. This is not the same as summing the series at its singularity.
* Periodic length-one and length-two conventions change tiny finite lattices but have no role in an
  increasing-side-length thermodynamic limit. They still matter critically for cross-checks and
  prefactors, which is why the discovered `c=1` bug was material.

### 10. Pauli-centraliser scope attack: theorem survives, old consequence fails

The elementary theorem inside the Pauli group is correct. But linear combinations of Pauli strings
can have commutator cancellations, and non-Pauli graph symmetries are not covered. The dense `C3`
translation certificate has

```
max |[U,A]| = max |[U,B]| = max |U^3-I| = 0,
det Gram(I, prod X, U) = 480 != 0.
```

This forced the repository to delete “only two invariant sectors” and state the actual conclusion:
there is no extensive family of commuting same-Hilbert-space Pauli-string fluxes. Momentum sectors
remain, and non-Pauli/ancilla mechanisms remain open.

### 11. Kac–Ward quantifier attack: exact subspace result survives, broad wording fails

The planar control is strong: modular arithmetic in `Z[zeta_8]`, CRT beyond a Hadamard bound and
direct even-subgraph enumeration agree coefficient-for-coefficient on four free square boxes.
For the strict cubic-covariant 3D scalar slice, straight and orthogonal transitions have only
weights `(a,b)`. The exact trace equations at lengths 4, 6 and 8 have Gröbner basis `{1}` only
after the order-eight equation is added. This genuinely rules out that two-parameter slice.

It does **not** rule out the full translation-invariant scalar family. That family has 30 weights,
25 after gauge fixing on one nonzero chart; its equations are consistent through order six and
the order-eight Gröbner computation timed out.

The audited text originally said:

```text
in 3D no such relationship holds for any scalar turn rule
```

and headed the comparison simply “3D (the relation cannot hold).” Both quantifiers were invalid.
After the audit, the sentence reads:

```text
in 3D no such relationship holds for any strict cubic-covariant scalar turn rule
```

and the heading carries the same qualifier, with an explicit pointer to the unresolved 30-weight
family. The numerical unsigned-half-angle failure at `v^6` is only a falsification of that natural
guess; it is not a general determinant no-go.

### 12. Cross-theorem synthesis attack: one obstruction had been conflated with another

The final report originally said all four no-go results “all point at the same object,” the local
quartic `Y_v Z_u Z_u' Z_u''`. This cannot include Theorem 11: its two-line Pauli-centraliser proof
does not use degree, claws, or a quartic term and holds unchanged on a 1D cycle, where the DG
quartic density is exactly zero. The corrected report groups only the claw/Gaussian, DG and TD
obstructions around `YZZZ`; it identifies the Pauli-flux result as a logically distinct no-go.

### 13. Artifact/provenance attack: current compliance failure

The current result contract requires `{"provenance":...,"data":...,"checks":[...]}`. Ten legacy
claim-bearing JSON files do not satisfy it:

```
results/dg_local_structure.json
results/dg_theorem_check.json
results/dla_table.json
results/dolan_grady.json
results/minimal_obstruction.json
results/onsager_algebra_dims.json
results/tridiagonal.json
results/tridiagonal_certificate.json
results/tridiagonal_certificate_general.json
results/reference_data/claimed_solutions.json
```

Most are raw lists; `dg_local_structure.json` is an unwrapped dictionary. This does not make their
integer contents false, but it weakens reproducibility and violates the repository's own stable API.
It also makes a blanket “all checks passed” count ambiguous because several artifacts have no check
objects at all.

## What survived and what did not

### Survived the attacks performed

* Theorem DG's closed formula, exactly, by dense integer matrices.
* The no-TD2 result on the `3x3` torus, with a large exact minimum residual.
* Both FLM geometric order bounds, after extending the arguments to Mayer clusters.
* Six low-order 3D series coefficients by direct non-FLM torus counting.
* Current periodic length-1/length-2 bond conventions and the repaired plus-boundary path.
* Safe directions and decimal containment of the certified `K_c` interval.
* Finite-volume duality's distinction between spin and gauge theories.
* The Pauli-group centraliser theorem at its corrected scope.
* The planar Kac–Ward control and the strict two-parameter cubic-covariant obstruction.

### Failed or required correction

* The original DG “iff 2-regular” corollary.
* The claim that every occurrence of the 3D benchmark was evaluation-only: the old Lee–Yang fit
  selection was benchmark-dependent.
* The old `plus_boundary=True, c=1` transfer polynomial.
* The inference from a two-element Pauli centraliser to a two-element full commutant / only two
  invariant sectors.
* The old Kac–Ward sentence claiming the strict-slice result for “any scalar turn rule.”
* The old synthesis saying the degree-independent Pauli-flux theorem points to the quartic defect.
* Repository-wide JSON result-contract compliance.

### Not broken, but not promoted

* TD for untested sizes remains open.
* The full 30-weight Kac–Ward family, zero-weight gauge branches, projective/spin-structure rules
  and bond dimension at least two remain open.
* Exponential algebra growth remains conjectural.
* Series coefficients above the independently checked ranges remain machine results, not new
  theorems.
* CFT numbers and several rigorous inequalities remain external inputs.
* Numerical roots, spectra and approximants remain numerical even at 100 digits.

## Residual weakest points and concrete actions

1. **Higher series coefficients have no fully independent implementation.** Extend direct GF(2)
   connected-subgraph enumeration from `v^10` to at least `v^12`, and direct droplet/Mayer
   enumeration from `x^12` to `x^16`; compare exact Fractions coefficient by coefficient.
2. **TD has a finite-size scope.** Derive a translation-local witness valid symbolically for every
   `Lx,Ly` above an explicit threshold, then treat collision cases (`L=2,3`) separately. Do not
   extrapolate the `3x3..5x5` table.
3. **The full commutant is unclassified.** Compute exact symmetry-reduced commutants for small
   layers and state separately what comes from space-group symmetry, global parity and accidental
   degeneracy. Theorem 11 must remain labelled “inside the Pauli group.”
4. **The infrared endpoint rests on external normalisation-sensitive theorems.** Quote the exact
   theorem statement and Fourier convention from the primary source beside each factor of two;
   independently enclose the Watson integral from its defining positive integral, not only its
   gamma evaluation.
5. **Simon–Lieb long-box arithmetic is not purely exact.** Re-evaluate the best useful box using
   rationals at a rational `v` bracketing the root or an outward-rounded interval library, so the
   certificate does not depend on a hand-counted floating-operation bound.
6. **The duality code does not test periodic side length one.** Either prove the cubical-chain code
   intentionally excludes that degenerate complex or add explicit tests reconciling dropped Ising
   self-loops with gauge cells. Continue to test side length two with parallel links.
7. **The Kaufman sector convention was inferred against the answer.** Supply a self-contained
   Pfaffian/spinor derivation of the finite-torus signs, or obtain and quote a primary source with
   the exact convention. Numerical agreement alone does not make the sign derivation independent.
8. **Legacy JSON lacks provenance/check envelopes.** Regenerate the ten listed files through their
   scripts under the required schema, including interpreter, arithmetic, and explicit pass/fail
   checks. Do not hand-wrap stale data.
9. **External literature data are not local proofs.** Preserve the E tag on bootstrap exponents,
   amplitude ratios, Lee–Yang, GKS/Simon–Lieb/infrared inequalities and NP-completeness scope. Add
   page/equation-level citations where currently only a paper-level reference is given.
10. **Benchmark-blind analysis needs process evidence, not only code flow.** For future series or
    edge fits, commit the selection rule and input coefficient hash before revealing benchmark
    comparisons; emit both hashes in the result artifact.
11. **The general Kac–Ward family is unresolved.** Cover zero-weight gauge charts and continue the
    exact order-eight system with elimination or certified local inconsistency witnesses. Keep the
    strict two-weight Gröbner result separate from all 30-weight, projective, multi-determinant and
    auxiliary-matrix constructions.

## Reproduction

From the repository root:

```text
.venv/bin/python tests/test_audit_replication.py
```

Observed dense outputs:

```text
DG star K1,3: PASS; max deviation=0
DG open 2x3 grid: PASS; max deviation=0
TD 3x3 torus no-solution: PASS; minimum Frobenius residual=6129.21297495;
  relative=0.244665468426; exact RSS=488374272/13
```

The script also writes the exact parameters, residuals, direct series counts, boundary cases,
literal scan, centraliser counterexample and endpoint calculations to
`results/audit/audit_replication.json`.
