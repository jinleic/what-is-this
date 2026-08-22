# Calibrated Lee–Yang edge scaling

## Scope and external control

- **[EXTERNAL]** Cardy’s two-dimensional conformal analysis identifies the Yang–Lee edge with the non-unitary minimal model conventionally denoted \(M(2,5)\) and gives the density exponent \(\sigma=-1/6\).  This value is a control target only; it is not an input to any fit, size/rank-window choice, or optimizer initialization.  Primary source: J. L. Cardy, “Conformal Invariance and the Yang-Lee Edge Singularity in Two Dimensions,” *Physical Review Letters* **54**, 1354–1356 (1985), doi:[10.1103/PhysRevLett.54.1354](https://doi.org/10.1103/PhysRevLett.54.1354), abstract.
- **[COMPUTATION]** The calibration coupling is fixed in advance at the exact rational bond fugacity \(x=e^{-2K}=2/3\).  This is also one of the frozen three-dimensional rational-grid points, so a successful control would permit an unchanged application without interpolating in coupling.
- **[EXTERNAL]** The exact square-lattice transition is \(K_c=\tfrac12\log(1+\sqrt2)=0.4406867935\ldots\) (`onsager1944` in the source manifest).  **[COMPUTATION]** The fixed calibration value \(K=-\tfrac12\log(2/3)=0.2027325541\ldots<K_c\), so this is a genuine high-temperature Yang–Lee edge control.
- **[COMPUTATION]** Calibration uses open square boxes, matching the open-cube boundary family available through \(4^3\).  For each \(L\), the field polynomial has exact integer coefficients after multiplication by the common factor \(3^{2L(L-1)}\).

## Estimator and stability criterion — fixed before inspecting any 3D scaling fit

- **[COMPUTATION]** Write the positive zero angles as \(0<\theta_1(L)<\theta_2(L)<\cdots\), put \(N=L^2\), and fit the first six ranks to the integrated-density model
  \[
  \theta_j(L)=\theta_{\rm edge}+A\left(\frac{j-1/2}{N}\right)^p,
  \qquad \widehat\sigma=p^{-1}-1.
  \]
  The nonlinear least-squares residual is unweighted in angle; the fitted bounds are \(0\leq\theta_{\rm edge}\leq\min\theta_j\), \(A>0\), and \(1/2\leq p\leq2\).  No literature exponent is supplied to the optimizer.
- **[COMPUTATION]** The predeclared calibration window is the eight largest consecutive exactly reachable sizes.  With the current public CRT transfer guard this is \(L=8,9,\ldots,15\); \(L=16\) would require \(2^{16}(16^2+1)=16{,}842{,}752\) final-layer residues and exceeds the guard of \(10{,}000{,}000\).
- **[COMPUTATION]** The control passes only if all of the following hold: (i) the full-window estimate obeys \(|\widehat\sigma+1/6|\leq0.03\); (ii) each of the eight leave-one-size-out refits obeys the same \(0.03\) accuracy tolerance; (iii) the range of those eight refits has width at most \(0.03\); and (iv) every fit converges away from the imposed parameter bounds.  These thresholds and the size/rank window are not changed after comparison with \(-1/6\).
- **[COMPUTATION]** If any control condition fails, no three-dimensional edge exponent is fitted or reported.  Three-dimensional first-zero angles and finite-volume near-edge densities remain independent exact-polynomial data products.


## Exact 2D polynomial reach

- **[COMPUTATION]** `e34_lee_yang_scaling.py` computed the exact integer field polynomial at \(x=2/3\) for every open square \(L=4,\ldots,15\) by modular transfer and CRT.  The largest polynomial is the degree-225 \(15\times15\) polynomial: it has 226 coefficients, common clearing factor \(3^{420}\), largest-coefficient length 234 decimal digits, and coefficient-array SHA-256 `aded33e7343315261b254ebfc49f77dfd3bd71733fb936f3e8d9b52d3b4551aa`.  Its full coefficient array is stored in `results/lee_yang/scaling.json`.
- **[COMPUTATION]** Every reduced Chebyshev polynomial was isolated over the rationals with all roots simple and inside \([-1,1]\), then refined at 100 decimal digits plus 30 guard digits.  For \(L=15\), exact isolation finds all 112 reduced roots and the exact factor \(z+1\) supplies the 113th positive-arc zero.  Across the calibration sequence the largest coefficient-max-normalized residual at the refined unit-circle roots is \(9.54\times10^{-109}\).
- **[COMPUTATION]** The exact-polynomial first-zero sequence at the predeclared coupling is:

| \(L\) | \(N=L^2\) | \(\theta_1(L)\) |
|---:|---:|---:|
| 4 | 16 | 0.898966462828489 |
| 5 | 25 | 0.776846720381498 |
| 6 | 36 | 0.702899187173097 |
| 7 | 49 | 0.654670235644043 |
| 8 | 64 | 0.621491296051092 |
| 9 | 81 | 0.597725592991487 |
| 10 | 100 | 0.580152872863183 |
| 11 | 121 | 0.566819984599356 |
| 12 | 144 | 0.556484915427976 |
| 13 | 169 | 0.548326955233039 |
| 14 | 196 | 0.541786238644718 |
| 15 | 225 | 0.536470427810958 |

## 2D calibration outcome

- **[COMPUTATION]** The predeclared \(L=8,\ldots,15\), \(j=1,\ldots,6\) fit gives \(\widehat\sigma=+0.0730115417808\), not \(-1/6\).  Its absolute control error is \(0.239678208447\), almost eight times the allowed \(0.03\).
- **[COMPUTATION]** The eight leave-one-size-out estimates range from \(+0.0590793219939\) to \(+0.0849621312842\).  Their span, \(0.0258828092903\), passes the \(0.03\) stability-width test, and all nine fits converge away from their bounds; nevertheless every estimate fails the required accuracy test.  Stability around the wrong value is not calibration.
- **[COMPUTATION]** Consecutive eight-size windows ending at \(L_{\max}=11,12,13,14,15\) give central estimates \(0.19830,0.16138,0.12326,0.09392,0.07301\), respectively.  This is visible drift toward the external value but not its recovery at the reached sizes.
- **[COMPUTATION]** The first uncomputed size that could alter the decision is \(L=16\), but it is beyond the current exact-transfer guard.  No finite required size is certified by these data.  As a strictly conditional resource diagnostic, fitting the maximum leave-one-out control error over the five rolling windows to \(aL_{\max}^{-\omega}\) gives \(\omega=1.3540\) and projects the \(0.03\) accuracy threshold only near \(L_{\max}=72\).  The present \(L=15\) final-layer array has \(7{,}405{,}568\) entries; the same representation at \(L=72\) would have \(24{,}485{,}470{,}213{,}679{,}110{,}433{,}013{,}760\), about \(3.31\times10^{18}\) times as many.  Thus this conditional projection is not a “twice the compute” extension and is neither a guarantee nor evidence that the same correction law survives that far.
- **[COMPUTATION]** **The 2D control fails.**  Therefore no 3D edge exponent was fitted, inferred, quoted, or compared with a literature value.  In particular, `scaling.json` stores `null` for the 3D exponent.  This refusal is the mandatory calibration-gate outcome, not a root-precision failure.
- **[COMPUTATION]** The decided limitation is **the predeclared one-term estimator on the available \(L\leq15\) lattice sequence**; the exact roots themselves are not the limiting factor.  The calculation does not distinguish an asymptotically valid estimator obscured by large finite-size corrections from a structurally inadequate estimator.  It therefore does not justify saying that larger lattices alone would cure the failure, although the monotone rolling-window drift makes that possibility testable.

## Independent 3D exact-zero and density data

- **[COMPUTATION]** Independently of exponent extraction, the frozen exact open-cube coefficient data were re-evaluated or imported for \(L=2,3,4\) at all six rational fugacities.  Every reduced root was isolated exactly; the new angles agree with the frozen `e04` angles to better than \(10^{-68}\).  The finite-volume density is
  \[
  \rho_{j+1/2}(L)=\frac{1}{L^3[\theta_{j+1}(L)-\theta_j(L)]}.
  \]
  `scaling.json` records the first eight intervals for every available cube and coupling.
- **[COMPUTATION]** The exact-polynomial first-zero sequences and the first \(L=4\) density interval are:

| \(x\) | \(K=-\tfrac12\log x\) | \(\theta_1(2)\) | \(\theta_1(3)\) | \(\theta_1(4)\) | \(\rho_{3/2}(4)\) |
|---:|---:|---:|---:|---:|---:|
| \(7/10\) | 0.178337471969 | 1.178335815297 | 0.638283238770 | 0.422325413530 | 0.100966076512 |
| \(2/3\) | 0.202732554054 | 1.080152546703 | 0.535068907812 | 0.323616207214 | 0.100963699637 |
| \(9/14\) | 0.220916376140 | 1.014485152796 | 0.469091894421 | 0.263225750438 | 0.101668473765 |
| \(5/8\) | 0.235001814623 | 0.967559228814 | 0.423733309640 | 0.223399747536 | 0.102609315775 |
| \(3/5\) | 0.255412811883 | 0.905109348464 | 0.366094420005 | 0.175501665721 | 0.104599426636 |
| \(1/2\) | 0.346573590280 | 0.692138336662 | 0.203275835223 | 0.0719716886854 | 0.123398616229 |

- **[COMPUTATION]** At every sampled coupling, \(\theta_1\) strictly decreases from \(L=2\) through \(L=4\).  The table is a finite-size trend only: it does not establish a thermodynamic edge, its limiting angle, or an exponent.

## Certified high-temperature finite-volume check

- **[THEOREM]** For the exact finite field polynomial \(P_L(z)=\sum_k A_kz^k\), every \(A_k\) is a positive integer at the tested ferromagnetic rational couplings.  Hence \(P_L(1)=\sum_kA_k>0\) exactly, so \(\theta=0\) is not a zero; because the root set is finite, its first positive angle is strictly separated from zero.
- **[COMPUTATION]** This repository’s certified bound \(K_c\geq0.2074277114992039908436804465100887760027\) proves that the tested points \(x=7/10\) (\(K=0.178337471969\)) and \(x=2/3\) (\(K=0.202732554054\)) lie above \(T_c\).  At both points, exact \(P_L(1)>0\) and the isolated sequences above verify \(\theta_1(L)>0\) for \(L=2,3,4\).  This is only the requested finite-size shadow; no size-uniform positive lower bound is claimed.

## Reproduction

- **[COMPUTATION]** Generate the artifact with `.venv/bin/python experiments/e34_lee_yang_scaling.py`; validate it independently with `.venv/bin/python tests/test_lee_yang_scaling.py`.  The experiment prints a final `PASS` while explicitly recording that the scientific calibration condition is false.