# RG and CFT constraints on a candidate exact 3D Ising solution

## Scope and bottom line

Renormalisation-group (RG) and conformal-field-theory (CFT) data provide stringent **necessary conditions** on the critical singularities of a proposed simple-cubic Ising free energy. They do not provide the microscopic free energy, and they do not determine the lattice coupling at which criticality occurs.

The machine-readable result is `results/cft/cft_constraints.json`; it was produced by `experiments/e25_cft_constraints.py` with 50-decimal-digit `mpmath` arithmetic [7]. The numerical singularity extractor:

* finds the positive-real curvature peak without being given a claimed critical coupling;
* correctly classifies the exact two-dimensional Onsager control as a logarithmic specific-heat singularity, hence `alpha = 0` [4,7];
* classifies both evaluable claimed three-dimensional free energies in `results/reference_data/claimed_solutions.json` as logarithmic, incompatible with the conservative three-dimensional constraint `alpha = 0.110(1)` [2,5-7]; and
* records the Rosengren entry as not testable by this checker because that record supplies a critical-point conjecture but no free-energy function [5,7].

This is a constraint/falsification result. It is not an exact solution of the lattice model, and neither a bootstrap island nor agreement with the data below would by itself constitute one.

## 1. Universal data used as constraints

### 1.1 Uncertainty convention

Kos, Poland, Simmons-Duffin, and Vichi find the numerical-bootstrap island

\[
\Delta_\sigma=0.5181489(10),\qquad
\Delta_\epsilon=1.412625(10)
\]

and quote the translations `eta = 0.0362978(20)` and `nu = 0.629971(4)` [1, abstract, introduction, and Eqs. (3.1)-(3.2)]. Their parenthesized ranges are numerical allowed-region/error estimates, not repeated-sampling Gaussian standard deviations. The derived uncertainties below are therefore **conservative rectangular propagations** of the two quoted marginal ranges. No unreported covariance is invented [1,7].

The leading correction exponent is a separate input. The downloaded precision-island paper does not determine the leading irrelevant scalar. The locally archived primary source that does quote it gives `omega = 0.83(5)`, equivalently the confluent exponent `Delta = omega*nu = 0.52(3)` [2, Table 1]. Its much larger uncertainty is retained rather than silently replaced by an uncited newer value.

### 1.2 Recommended constraint set

| quantity | value used | status and provenance |
|---|---:|---|
| `Delta_sigma` | `0.5181489(10)` | Numerical-bootstrap island [1, Eq. (3.1)]. |
| `Delta_epsilon` | `1.412625(10)` | Numerical-bootstrap island [1, Eq. (3.2)]. |
| `nu` | `0.629971(4)` | `1/(3-Delta_epsilon)`; the same rounded value is quoted by Kos et al. [1, introduction; 7]. |
| `eta` | `0.0362978(20)` | `2 Delta_sigma-1`; also quoted by Kos et al. [1, introduction; 7]. |
| `omega` | `0.83(5)` | Leading correction-to-scaling exponent from the standard-model HT-series analysis [2, Table 1]. |
| `alpha` | `0.110087(12)` | `2-3 nu`, with rectangular propagation of the bootstrap range [1,7]. For the checker we deliberately use the wider `0.110(1)` constraint requested in this track and supported by the direct/derived series values in [2, abstract and Table 1]. |
| `beta` | `0.326419(3)` | `nu(1+eta)/2`, with rectangular propagation [1,7]. |
| `gamma` | `1.237075(9)` | `nu(2-eta)`, with rectangular propagation [1,7]. |
| `delta` | `4.789841(12)` | `(5-eta)/(1+eta)`, with rectangular propagation [1,7]. |

These last four high-precision values are scaling consequences of the two bootstrap dimensions. Calling them independent measurements would be wrong.

## 2. Scaling-relation consistency checks

For `d = 3`, the checked identities are

\[
\alpha=2-3\nu,\qquad
\gamma=\nu(2-\eta),\qquad
\beta={\nu(1+\eta)\over2},\qquad
\delta={5-\eta\over1+\eta}.
\]

Substitution of the unrounded values derived from `(Delta_sigma, Delta_epsilon)` gives a largest stored residual of `4.44e-32`, the expected decimal-serialization roundoff [1,7]. That check catches transcription and formula errors, but it is algebraic rather than statistically independent.

A more useful cross-source test compares those predictions with separately reported thermodynamic estimates. The “pull” below is only a diagnostic ratio formed by adding displayed uncertainties in quadrature: the bootstrap ranges are not Gaussian errors, covariances are unavailable, and the final `delta` source explicitly marks its value as scaling-derived [1-3,7].

| relation | published left side | bootstrap-dimension prediction | residual `lhs-rhs` | diagnostic absolute pull | sources |
|---|---:|---:|---:|---:|---|
| `alpha = 2-3 nu` | `0.110(2)` | `0.1100874085(119)` | `-8.74085e-5` | `0.0437` | Direct HT `alpha` in [2, Table 1]; prediction [1,7]. |
| `gamma = nu(2-eta)` | `1.23708(33)` | `1.237075171(91)` | `+4.82873e-6` | `0.0146` | Monte Carlo [3, Eq. (46)]; prediction [1,7]. |
| `beta = nu(1+eta)/2` | `0.32630(22)` | `0.326418710(27)` | `-1.18710e-4` | `0.540` | Monte Carlo [3, Eq. (47)]; prediction [1,7]. |
| `delta = (5-eta)/(1+eta)` | `4.7893(8)` | `4.78984149(12)` | `-5.41492e-4` | `0.677` | Reported scaling-derived value [2, Table 1]; prediction [1,7]. |

No inconsistency is seen: every central-value residual is smaller than one displayed combined uncertainty. The largest diagnostic ratio is `0.677` [7]. This does not make the inputs exact; it shows that the independently obtained Monte Carlo/series numbers do not conflict with the bootstrap scaling data at their quoted precision.

## 3. What the data constrain—and what they do not

### They constrain

1. **The leading singular part.** With `t=(T-T_c)/T_c`, a three-dimensional Ising free energy must have a singular contribution consistent with
   \[
   f_s(t,0)\sim F_\pm |t|^{2-\alpha},\qquad \alpha=0.110087(12),
   \]
   up to corrections governed first by `omega = 0.83(5)` and analytic terms [1,2]. For this repository's reduced free energy `phi=(1/N) log Z`, the dimensionless specific heat is proportional to `K^2 phi''(K)`. The analytic factor `K^2` does not change the exponent at nonzero `K_c`.
2. **Universal critical exponents and scaling dimensions.** Correlation length, spin correlations, magnetisation, susceptibility, and the critical isotherm must reproduce the values in the table above [1,2].
3. **Universal ratios/combinations of amplitudes.** Individual metric factors and amplitudes depend on microscopic normalisations, but ratios defined with consistent thermal and magnetic conventions are universal [2, Sec. 5].
4. **The exponent of corrections, not their coefficient.** `omega` is universal; the amplitude multiplying the leading irrelevant correction is not and can vanish in an improved lattice Hamiltonian [2, Secs. 1 and 3].

### They do not constrain

1. **`K_c` is non-universal.** A CFT describes the fixed point, not the non-universal map from the simple-cubic lattice coupling `K=beta J` to the thermal scaling field. The Monte Carlo value `K_c=0.221654626(5)` is an independent lattice benchmark [3, abstract], not a CFT prediction.
2. **The regular background is not fixed.** Adding a function analytic through `K_c` changes `phi`, energy, and the background part of the specific heat without changing the CFT singular data.
3. **Absolute thermodynamic amplitudes are generally not fixed.** Rescaling the thermal or magnetic scaling fields changes them. Only specified universal ratios or combinations survive [2, Sec. 5].
4. **A zero-field free energy cannot by itself expose all magnetic data.** Testing `beta`, `gamma`, `delta`, `Delta_sigma`, `C_+/C_-`, or `R_chi` requires field derivatives, spontaneous magnetisation, or correlation functions, not merely `phi(K,0)`.
5. **A numerical bootstrap island is not an exact lattice solution.** The island is a numerically isolated set of CFT data under stated spectrum assumptions and finite numerical searches [1, Secs. 2-3]. It does not produce finite-volume partition functions, the lattice regular background, the map from `K` to the scaling field, or `K_c`.

## 4. The constraint checker

The public function

```python
check_candidate(phi_callable, name)
```

accepts a callable for the reduced zero-field free energy `phi(K)`. It does not accept or use a proposed `K_c`. Its recorded procedure is [7]:

1. scan the positive-real interval `0.03 <= K <= 0.58` for significant peaks of the centered second difference of `phi`;
2. choose the nearest significant peak and refine its position on seven successively finer grids;
3. evaluate a five-point estimate of `phi''` on both sides at reduced offsets from `0.08` down to `0.000625`;
4. remove an analytic constant background by differencing successive curvatures;
5. use the ratio of those differences under offset halving. For `C_s ~ |t|^-alpha`, the ratio tends to `2^alpha`; for `C_s ~ log|t|`, the differences tend to a nonzero constant and the effective exponent tends to zero; and
6. quote an error equal to the largest of the low/high-side asymmetry, last-level drift, singularity-location contribution, and a numerical floor. It also compares a logarithmic fit with a fixed-`alpha=0.110` fit [2,7].

This is a numerical diagnostic over the stated interval, not a proof that no other complex-plane singularity exists. Its location estimate is deliberately not used as a universality constraint, because `K_c` is non-universal.

### 4.1 Onsager control

For the exact isotropic square-lattice free energy implemented in `ising.onsager.onsager_free_energy`, the checker finds

* `K_s = 0.4406873 +/- 0.0000024`, whose interval contains the exact `K_c = log(1+sqrt(2))/2 = 0.4406867935097715...` [4,7];
* raw two-sided effective `alpha = 0.00024 +/- 0.01034`; and
* classification **`logarithmic_alpha_0`** [7].

The mean relative RMS is `0.00375` for the logarithmic model versus `0.01863` for a fixed `alpha=0.110` power [7]. Thus the extractor detects the known logarithmic singularity rather than confusing a small positive exponent with zero. A synthetic `alpha=0.110` power-law control is separately classified as a power law, with extracted `alpha = 0.110001 +/- 0.01644` [2,7].

### 4.2 Claimed three-dimensional formulas

The critical couplings shown in the next table were **not** supplied to the locator. Their agreement with each claim merely verifies that the curvature peak being analysed is the intended one; it is not a CFT test [5-7].

| record | nearest detected positive-real singularity | extracted specific-heat behavior | fit comparison: log vs fixed `alpha=0.110` | verdict |
|---|---:|---:|---:|---|
| Zhang 2007 | `0.2406056 +/- 0.0000024` [5,7] | raw `alpha = 0.00030 +/- 0.01595`; **logarithmic** [7] | `0.00429` vs `0.01856` [7] | Incompatible with `alpha=0.110(1)` [2,7]. |
| Degang Zhang 2021 | `0.3046895 +/- 0.0000024` [6,7] | raw `alpha = 0.00030 +/- 0.01259`; **logarithmic** [7] | `0.00120` vs `0.01844` [7] | Incompatible with `alpha=0.110(1)` [2,7]. |
| Rosengren critical-point conjecture | Claimed `0.221658637208699`, but no `phi(K)` is stored [5,7] | Not extractable | Not applicable | `not_testable_no_free_energy_callable`; no singular exponent is inferred [7]. |

For the Zhang 2007 expression, the script evaluates one of the two angular integrals analytically using the standard integral of `log(A-B cos theta)`; this is an exact dimensional reduction of the expression stored in `claimed_solutions.json`, not a fitted surrogate [5,7]. The 2021 expression is transcribed directly with cancellation-free evaluation near its branch point [6,7]. Source-expression hashes and all curvature sequences are stored in the result JSON [7].

These failures are decisive for the recorded formulas as candidate simple-cubic free energies: a logarithmic `alpha=0` singularity is outside even the deliberately widened three-dimensional interval `0.110(1)` [2,7]. They do not prove that no different exact representation can exist.

## 5. Universal amplitude ratios

Campostrini et al. define the specific-heat amplitudes by `C_H=A_+ t^-alpha` above and `C_H=A_- |t|^-alpha` below the transition, the susceptibility amplitudes by `chi=C_+ t^-gamma` and `chi=C_- |t|^-gamma`, and the magnetisation amplitudes `B` and `B_c` on the coexistence curve and critical isotherm [2, Eqs. (4.3), (4.4), (4.10), (5.1), and Table 7]. With those conventions their estimates are:

| ratio | value | definition and source |
|---|---:|---|
| `A_+/A_-` | `0.532(3)` | `U0`; [2, Tables 1 and 7]. |
| `C_+/C_-` | `4.76(2)` | `U2`; [2, Tables 1 and 7]. Here `C_+/-` denote susceptibility amplitudes, not heat capacities. |
| `R_chi` | `1.660(4)` | `C_+ B^(delta-1)/B_c^delta`; [2, Tables 1, 7, and 8]. |

These are sourced numerical estimates, not exact constants and not outputs of the bootstrap island [2]. A future candidate should reproduce them only after matching the same definitions of `t`, field, magnetisation, and susceptibilities.

### Can the repository's exact series test them now?

No. The current extended artifacts (`results/series/extended_sc_ht_free_energy.json`,
`results/series/extended_sc_lt_free_energy.json`, produced by `experiments/e18_series_extend.py`)
contain the zero-field HT free energy through `v^20` -- **ten** nonzero coefficients
(`v^2, v^4, ..., v^20`) -- and the zero-field LT free energy through `x^28` -- **eleven** nonzero
coefficients (`x^6, x^10, x^12, ..., x^28`). The earlier baseline artifacts from
`experiments/e03_flm_series.py` stopped at `v^16` and `x^20`; the counts below were updated to the
extended artifacts after an adversarial audit flagged the stale inventory. The conclusion is
unchanged. Specifically:

* `A_+/A_-` would require controlled analytic continuation to the same critical point from both
  phases and separation of the weak `alpha approximately 0.11` singularity from regular
  backgrounds. Ten nonzero HT terms and eleven nonzero LT terms do not support an honest amplitude
  estimate [7,8]; for comparison, `experiments/e19_series_analysis.py` extracts only
  `alpha = 0.35 +- 0.34` from the same HT series, an uncertainty three times the value.
* `C_+/C_-` requires susceptibility series on both sides; those field derivatives are absent [2,7,8].
* `R_chi` additionally requires `B` and `B_c`, so it is also absent [2,7,8].

Therefore the set of amplitude ratios testable by the present exact series is **none** [7,8]. Inferring any of them from these short zero-field series would be overfitting.

## 6. Concrete checklist for a future candidate

A proposed exact simple-cubic Ising solution should pass all applicable items below. Passing them is necessary, not sufficient.

1. **Lattice identity and normalisation.** It must compute the model with `K=beta J` and the partition-function conventions in `problem_specification.md`, including boundary/topological sectors before taking the thermodynamic limit.
2. **Independent non-universal location.** Its nearest physical positive-real singularity should reproduce the lattice benchmark `K_c=0.221654626(5)` [3, abstract]. This is a lattice test, explicitly not a CFT prediction.
3. **Thermal singularity.** `K^2 phi''(K)` must have a power-law singular part with `alpha=0.110087(12)` [1,7], or at minimum overlap the conservative test interval `0.110(1)` [2,7]. A logarithm (`alpha=0`) fails.
4. **CFT operator dimensions.** It must be consistent with `Delta_sigma=0.5181489(10)` and `Delta_epsilon=1.412625(10)` [1].
5. **Correlation data.** It must give `nu=0.629971(4)` and `eta=0.0362978(20)` [1,7].
6. **Magnetic singularities.** With a field or equivalent observables, it must give `beta=0.326419(3)`, `gamma=1.237075(9)`, and `delta=4.789841(12)` [1,7].
7. **Corrections to scaling.** Generic leading corrections must use `omega=0.83(5)`; their amplitude remains model-dependent and may vanish in an improved model [2, Table 1].
8. **Scaling identities.** The four identities in Section 2 must hold before rounded display. Independent observables should remain consistent with them; the present cross-source diagnostic ratios are all below `0.677` [1-3,7].
9. **Amplitude ratios.** With the conventions of [2], it should reproduce `A_+/A_-=0.532(3)`, `C_+/C_-=4.76(2)`, and `R_chi=1.660(4)` [2]. A zero-field formula can address only the first, and only if both phases and backgrounds are controlled.
10. **Regular and finite-volume information.** It must reproduce exact HT/LT coefficients and finite-volume partition functions as separate lattice checks. Universal critical data cannot substitute for these non-universal requirements.
11. **No circular fit.** None of the values above may be inserted as a fitted parameter and then presented as a prediction.

## Reproduction

From the repository root:

```text
.venv/bin/python experiments/e25_cft_constraints.py
```

The command prints the individual checks followed by `PASS` and rewrites `results/cft/cft_constraints.json` [7].

## References

1. F. Kos, D. Poland, D. Simmons-Duffin, and A. Vichi, *Precision Islands in the Ising and O(N) Models*, JHEP **08** (2016) 036, [doi:10.1007/JHEP08(2016)036](https://doi.org/10.1007/JHEP08(2016)036), [arXiv:1603.04436](https://arxiv.org/abs/1603.04436). Local full text: `sources/fulltext/num_kos2016.pdf`.
2. M. Campostrini, A. Pelissetto, P. Rossi, and E. Vicari, *25th-Order High-Temperature Expansion Results for Three-Dimensional Ising-Like Systems on the Simple-Cubic Lattice*, Phys. Rev. E **65** (2002) 066127, [doi:10.1103/PhysRevE.65.066127](https://doi.org/10.1103/PhysRevE.65.066127), [arXiv:cond-mat/0201180](https://arxiv.org/abs/cond-mat/0201180). Local full text: `sources/fulltext/num_campostrini2002.pdf`.
3. A. M. Ferrenberg, J. Xu, and D. P. Landau, *Pushing the Limits of Monte Carlo Simulations for the Three-Dimensional Ising Model*, Phys. Rev. E **97** (2018) 043301, [doi:10.1103/PhysRevE.97.043301](https://doi.org/10.1103/PhysRevE.97.043301), [arXiv:1806.03558](https://arxiv.org/abs/1806.03558). Local full text: `sources/fulltext/num_ferrenberg2018.pdf`.
4. L. Onsager, *Crystal Statistics. I. A Two-Dimensional Model with an Order-Disorder Transition*, Phys. Rev. **65** (1944) 117-149, [doi:10.1103/PhysRev.65.117](https://doi.org/10.1103/PhysRev.65.117). Metadata: `sources/manifest.yaml`; local independent implementation: `src/ising/onsager/`.
5. Z.-D. Zhang, *Conjectures on Exact Solution of Three-Dimensional Simple Orthorhombic Ising Lattices*, Philosophical Magazine **87** (2007) 5309-5419, [doi:10.1080/14786430701646325](https://doi.org/10.1080/14786430701646325), [arXiv:0705.1045](https://arxiv.org/abs/0705.1045). Local full text: `sources/fulltext/cla_zhang2007.pdf`.
6. D. Zhang, *Exact Solution for Three-Dimensional Ising Model*, Symmetry **13** (2021) 1837, [doi:10.3390/sym13101837](https://doi.org/10.3390/sym13101837), [arXiv:2110.11233](https://arxiv.org/abs/2110.11233). Local full text: `sources/fulltext/cla_zhang2021.pdf`.
7. Machine result generated here: `results/cft/cft_constraints.json`, provenance `experiments/e25_cft_constraints.py`, 50-decimal-digit arithmetic.
8. Exact local series artifacts: `results/series/sc_ht_free_energy.json` (`v^16`) and `results/series/sc_lt_free_energy.json` (`x^20`), generated by `experiments/e03_flm_series.py` with exact integer/rational arithmetic.
