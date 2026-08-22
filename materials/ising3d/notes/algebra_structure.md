# Abstract structure of the two-sum Ising Lie algebras

## 1. Scope and arithmetic

For a finite layer graph \(\Gamma=(V,E)\), put
\[
A=\sum_{v\in V}X_v,\qquad B=\sum_{(u,v)\in E}Z_uZ_v,
\qquad \mathfrak g_\Gamma=\operatorname{Lie}\langle A,B\rangle.
\]

[COMPUTATION] `experiments/e29_algebra_structure.py` builds every generated basis with
`ising.clifford.fast_lie.GeneratedAlgebra` in Pauli-string coordinates.  The \(2\times2\)
calculation uses exact `Fraction` arithmetic over \(\mathbb Q\).  Larger calculations use exact
arithmetic over
\(p_1=2147483647\) and \(p_2=2147483629\).  A nonzero modular minor is a rigorous lower bound for
the corresponding rational rank; agreement at two primes is a cross-check, not a
characteristic-zero equality proof.  Every larger-grid structural statement below is therefore
explicitly scoped to \(\mathbb F_p\).

[COMPUTATION] The machine records are `results/algebra_structure/structure.json` and
`results/algebra_structure/onsager_controls.json`.  Fields produced in the run name their method.
Characteristic-zero Killing forms, solvable radicals, and Levi decompositions are outside this
experiment's scope and are recorded as `NOT_COMPUTED_BY_THIS_EXPERIMENT`, never inferred from
modular image data.  Each grid record also carries a `characteristic_zero_resolution` pointer
naming the artifact (if any) that owns the exact resolution: for the \(2\times3\) layer that
owner is `results/algebra_structure/char0_levi.json` (`experiments/e45_char0_levi.py`); for the
\(2\times4\) layer no such artifact exists yet.

## 2. Results at a glance

| layer | generated dimension | computed centre | computed derived rank | abstract structure |
|---|---:|---:|---:|---|
| \(2\times2=C_4\) | 11 exactly over \(\mathbb Q\) | 2 exactly | 9 exactly; the derived ideal is \(A_1^3\) | [COMPUTATION] \(A_1^3\oplus\mathbb C^2\) |
| open \(2\times3\) | 263 at both primes; hence \(\dim_\mathbb Q\mathfrak g\ge263\) | 1 at both primes | 262 at both primes; hence rational rank \(\ge262\) | [COMPUTATION] \(\mathbb F_p\) only here; resolved exactly over \(\mathbb Q\) by e45: radical = centre (dim 1), Levi \(C_7\oplus C_3\oplus C_3\oplus A_8\oplus A_5\) (`results/algebra_structure/char0_levi.json`, `proofs/char0_structure.md`) |
| open \(2\times4\) | 2952 at both primes; hence \(\dim_\mathbb Q\mathfrak g\ge2952\) | 1 at both primes | 2951 at both primes; hence rational rank \(\ge2951\) | [COMPUTATION] characteristic-zero Killing form, radical, and Levi decomposition undetermined |

[COMPUTATION] For each modular case, the centre is computed as
\(\ker(\operatorname{ad}_A|\mathfrak g)\cap\ker(\operatorname{ad}_B|\mathfrak g)\).  This is the
centre because \(A,B\) generate \(\mathfrak g\).  The derived rank is computed as
\[
\dim\bigl(\operatorname{ad}_A(\mathfrak g)+
          \operatorname{ad}_B(\mathfrak g)\bigr)=\dim[\mathfrak g,\mathfrak g].
\]
Only the first two terms of the derived and lower-central series are recorded for the modular
grids.  No assertion that the derived algebra is perfect is made.

## 3. Exact-\(\mathbb Q\) four-cycle result

[COMPUTATION] The full exact rational \(11\times11\) Killing matrix is stored in
`structure.json`.  Its rank is 9 and its nullspace has dimension 2.  The simultaneous exact
kernels of \(\operatorname{ad}_A\) and \(\operatorname{ad}_B\) have dimension 2, and direct
commutators show that this space is central.  The exact derived ideal has dimension 9 and is
perfect.

[COMPUTATION] Resolving spin flip, row reflection, and column reflection gives sector-image
dimensions
\((6,4,0,1,0,1,0,4)\).  Joint image ranks distinguish three independent
three-dimensional simple ideals; each has Cartan matrix \([2]\).  Thus
\[
\boxed{\mathfrak g_{C_4}\cong
\mathfrak{sl}_2\oplus\mathfrak{sl}_2\oplus\mathfrak{sl}_2\oplus\mathbb C^2.}
\]
The dimension identity is \(3+3+3+2=11\).  Therefore the naive
\(\mathfrak{sp}_4\oplus\mathbb C\) or \(\mathfrak{so}_5\oplus\mathbb C\) guess is false: its
centre would have dimension one and its root system would be irreducible \(B_2=C_2\), whereas the
computed root system is \(A_1^3\).

## 4. Open \(2\times3\): finite-field quotient images only

[COMPUTATION] At both primes the generated dimension is 263, centre dimension is 1, and derived
rank is 262.  These are direct closure/kernel/rank computations over \(\mathbb F_p\), not claims
about the characteristic-zero radical or Killing form.

[COMPUTATION] The eight symmetry-sector module dimensions are
\((14,10,6,6,6,10,6,6)\), and their Lie-image dimensions are
\((105,81,21,36,21,81,21,36)\).  The following image identifications are certified separately at
each prime by saturation:

* the 105-dimensional image on 14 dimensions preserves a unique nondegenerate alternating form,
  so it equals \(\mathfrak{sp}_{14}(\mathbb F_p)\), type \(C_7\);
* each 21-dimensional image on 6 dimensions preserves a unique nondegenerate alternating form,
  so it equals \(\mathfrak{sp}_6(\mathbb F_p)\), type \(C_3\);
* each 36-dimensional image on 6 dimensions saturates \(\mathfrak{gl}_6(\mathbb F_p)\).

[COMPUTATION] The 81-dimensional image on a 10-dimensional module preserves neither a symmetric
nor an alternating form, but it does not saturate \(\mathfrak{gl}_{10}\).  Its Cartan type is
**not identified**.  Pairwise joint-image ranks do not determine kernels and their intersections,
so they do not yield a direct-sum decomposition of \(\mathfrak g\) **within this experiment**.
The characteristic-zero resolution was subsequently obtained by `experiments/e45_char0_levi.py`
(full-orbit exact-\(\mathbb Q\) closure, Killing rank 262, radical = centre of dimension 1, Levi
type \(C_7\oplus C_3\oplus C_3\oplus A_8\oplus A_5\)); see `proofs/char0_structure.md`.  The
\(\mathbb F_p\)-only statements of this section remain correct as statements about e29's own
computation.

## 5. Open \(2\times4\): partial finite-field data and wall

[COMPUTATION] At both primes the generated dimension is 2952, centre dimension is 1, and derived
rank is 2951.  The eight symmetry-sector module dimensions are
\((44,32,28,32,28,32,28,32)\), with image dimensions
\((861,606,351,799,351,799,351,606)\).

[COMPUTATION] On the `000` sector, a two-dimensional kernel leaves a 42-dimensional faithful
quotient with a unique nondegenerate symmetric form.  Its 861-dimensional image saturates
\(\mathfrak{so}_{42}(\mathbb F_p)\), type \(D_{21}\).  On the 28-dimensional sectors, a
one-dimensional kernel leaves a 27-dimensional faithful quotient with a unique nondegenerate
symmetric form.  Its 351-dimensional image saturates \(\mathfrak{so}_{27}(\mathbb F_p)\), type
\(B_{13}\).  These are explicitly \(\mathbb F_p\)-only quotient-image statements.

[COMPUTATION] The remaining 606- and 799-dimensional image types are unidentified.  All sector
maps are non-faithful and overlap; their visible image factors cannot be added as independent
ideals.  The abstract Killing form, solvable radical, Levi dimension, rank, and type are therefore
**UNDETERMINED**.

[COMPUTATION] The exact resource wall is the full adjoint/ideal decomposition.  Retaining 2952
dense \(2952\times2952\) int64 adjoint matrices would require about 205.8 GB.  The implemented
ambient-Pauli kernel and rank computations avoid that allocation and determine the centre and first
derived rank, but they do not recover the ideal-kernel intersection lattice needed for a Levi
decomposition.

## 6. One-dimensional Onsager controls

[THEOREM] The universal Onsager algebra is
\[
\mathrm{OA}=(\mathfrak{sl}_2\otimes\mathbb C[t,t^{-1}])^\theta,
\]
where \(\theta\) combines the Chevalley involution with \(t\leftrightarrow t^{-1}\).  The
Dolan--Grady relations give the homomorphism from OA to every chain/ring spin image.

[EXTERNAL] The exact finite quotient names used here are representation-theoretic:
\(\mathrm{OA}/\ker\rho_{\mathrm{open},n}\), where \(\rho_{\mathrm{open},n}\) is the length-\(n\)
Jordan--Wigner spin representation, and
\(\mathrm{OA}/\ker\rho_{\mathrm{per},n}\), where \(\rho_{\mathrm{per},n}\) is the periodic spin
representation including its fermion-parity boundary sectors.  This does not assert an unproved
principal-polynomial presentation of either kernel.

[COMPUTATION] For \(2\le n\le6\), the open-chain images have dimension \(n^2\), centre dimension
1, and derived rank \(n^2-1\) at both primes.  For \(3\le n\le6\), the ring images have dimension
\(3n-1\), centre dimension 2, and derived rank \(3n-3\) at both primes.  These reproduce the
expected Onsager spin-representation controls.  As for the grids, later derived-series terms are
not inferred from the first derived rank.

## 7. What the comparison does and does not establish

[COMPUTATION] The exact-\(\mathbb Q\) cycle control has the fully classified structure
\(A_1^3\oplus\mathbb C^2\).  Genuine two-dimensional layers have much larger modular closures and
large saturated classical quotient images over both tested finite fields.  This is a sharp
finite-field structural contrast, but it is not a characteristic-zero Levi classification.

[THEOREM] The \(2\times3\) layer algebra is **not** an abstract quotient of the Onsager
algebra: the exact characteristic-zero structure certified by e45 (radical = centre, Levi
\(C_7\oplus C_3\oplus C_3\oplus A_8\oplus A_5\)) combined with the closed-ideal Date--Roan
classification (the Levi factor of every finite-dimensional OA quotient is \(\mathfrak{sl}_2^n\))
excludes every realization; see `proofs/onsager_quotient_nogo.md` and `proofs/char0_structure.md`.

[CONJECTURE] For the \(2\times4\) layer the same question remains **OPEN**: no characteristic-zero
certification exists, and the modular \(D_{21},B_{13}\) image witnesses have no rational lift in
this experiment.  The repository's independent characteristic-zero obstruction to Onsager
integrability remains the proved Dolan--Grady defect theorem; it does not depend on this
exploratory classification computation.
