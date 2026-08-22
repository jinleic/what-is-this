# Graded interacting auxiliary \(R\) for the Ising tetrahedron relation

## Scope and reproduction

[COMPUTATION] This calculation addresses the unrestricted **grade-preserving** \(8\times8\)
auxiliary matrix in the fixed placement

\[
R_{123}L_{145}L_{246}L_{356}
 =L_{356}L_{246}L_{145}R_{123}.
\]

It does not assume that \(R\) is a direct sum, a tensor product, a gauge transform, or a product of
two-leg operators.  It also imposes no translation, cyclic-leg, or full leg-permutation symmetry:
every entry allowed by the grading is independent.  Thus genuinely three-leg interacting matrices
are included inside the stated graded subfamily.

[COMPUTATION] Reproduction commands are

```text
.venv/bin/python experiments/e44_tetra_graded.py
.venv/bin/python tests/test_tetra_graded.py
```

[COMPUTATION] The machine certificate is
`results/integrability/tetra_graded.json`.  Exact claims use Python integers and `Fraction`, and
SymPy polynomial arithmetic over \(\mathbb Q\).  Floating-point wall times are metadata only and do
not support any mathematical conclusion.

## Prerequisite: the physical six-leg tensor

[LEMMA] Put \(q=w^2=\tanh K\), and omit the uniform factor two at one vertex.  The verified Ising
site tensor, reshaped as an operator with three input and three output legs, is

\[
 L^{(q)}_{o,i}=\frac12T[i_1,i_2,i_3,o_1,o_2,o_3]
 =\begin{cases}
 q^{(|i|+|o|)/2},& |i|+|o|\equiv0\pmod2,\\
 0,& |i|+|o|\equiv1\pmod2.
 \end{cases}
\]

The omitted factor is harmless here: each side of RLLL contains exactly three copies of \(L\), so
restoring \(T=2L\) multiplies both sides by \(2^3\).

[COMPUTATION] `e44_tetra_graded.py` imports and reuses `finite_lattice_controls` and
`site_tensor_entry` from `experiments/e35_tetrahedron_spectral.py`.  Its artifact builder executes
these controls before constructing the RLLL system.  The exact tensor contraction again agrees
coefficient by coefficient with the independent finite-lattice enumerator and gives the required
scaled integers:

| open box | scaled partition integer |
|---|---:|
| \(2\times2\times2\) | 36,450 |
| \(2\times2\times3\) | 16,394,562 |
| \(2\times3\times3\) | 246,853,161,090 |

[COMPUTATION] At the exact interior specialization \(q=1/4\) (the Wave-5 choice \(w=1/2\)), all
64 local entries were also compared directly with the Wave-5 builder: there were zero mismatches.
The finite-box controls themselves retain their stored comparison point \(q=1/3\); their purpose is
to fix the tensor normalization and bond convention, not to tune the RLLL calculation.

[LEMMA] The physical operator preserves total \(\mathbb Z_2\) grade.  Indeed,
\(|i|+|o|\) is even exactly when \(|i|\equiv|o|\pmod2\).  The direct exact scan found 32 nonzero
entries, all in the 32 grade-allowed positions, and no violation.  In the even/odd ordered basis,
\(L\) has two \(4\times4\) diagonal grade blocks; at \(q=1/4\) each block has rank one and the full
matrix has rank two.

## Exact graded ansatz and polynomial system

[COMPUTATION] Let

\[
 {\cal C}=\{(o,i)\in\{0,1\}^3\times\{0,1\}^3:
 |o|\equiv|i|\pmod2\}.
\]

For every \((o,i)\in{\cal C}\), introduce an independent coordinate \(r_{o,i}\), and set every
other entry of \(R\) to zero.  There are

\[
 |{\cal C}|=4^2+4^2=32
\]

unknowns rather than 64.  Relative to total auxiliary parity, \(R\) is an arbitrary
\(4\times4\oplus4\times4\) matrix.  No relation is imposed among entries inside either block.
The equation is homogeneous of degree one in \(R\), so a nonzero solution has an overall scalar
gauge.  The primary normalization chart is \(r_{000,000}=1\); all 32 possible nonzero-coordinate
charts are checked below so this chart choice does not discard another projective solution.

[COMPUTATION] On the six binary spaces numbered 1 through 6, define

\[
 P(q)=L_{145}(q)L_{246}(q)L_{356}(q),\qquad
 Q(q)=L_{356}(q)L_{246}(q)L_{145}(q).
\]

The 4096 component equations are the entries of

\[
 F(r;q)=R_{123}P(q)-Q(q)R_{123}=0.
\]

Every component is a polynomial in \(q\) and a homogeneous linear form in the 32 coordinates of
\(R\).  At \(q=1/4\), the program uniformly replaces \(L\) by \(64L\), whose entries are integers
\(64,16,4,1\).  This merely multiplies every component equation by \(64^3\).

[COMPUTATION] Exact component extraction at \(q=1/4\) gives:

| quantity | exact value |
|---|---:|
| components | 4,096 |
| identically zero components | 2,048 |
| nonzero component linear forms | 2,048 |
| distinct primitive forms up to a rational scalar | 1,024 |
| independent graded unknowns | 32 |
| exact rational rank | 32 |
| nullity | 0 |

The exact row reduction selects 32 actual RLLL components.  Their complete integer coefficient
matrix and component labels are stored in the JSON artifact.  Its determinant is nonzero and
factors as

\[
 2^{297}3^{45}5^{39}.
\]

[COMPUTATION] Timings for the recorded run were 0.128527 s to construct the full component matrix
and 0.547431 s to extract and evaluate the exact rank minor.  As an independent implementation
check, the program embedded each of the 32 elementary local \(R\) matrices directly into the
\(64\times64\) six-space operator and compared all 131,072 coefficients with the derived component
matrix.  It found zero mismatches in 0.160180 s.

## Uniform parameter certificate

[LEMMA] Use the same 32 component rows before specializing \(q\).  Exact determinant evaluation in
\(\mathbb Q[q]\) gives

\[
 \det B(q)=-2q^{56}(q-1)^{45}(q+1)^{39}.
\]

Specializing this formula to \(q=1/4\) and multiplying each of the 32 rows by \(64^3\) reproduces
the fixed integer determinant above exactly.  The symbolic determinant computation took 1.316319
s in the recorded run.

[THEOREM] Let \(k\) be a characteristic-zero field and let
\(q\in k\setminus\{0,1,-1\}\).  The only \(\mathbb Z_2\)-grade-preserving \(8\times8\) matrix
satisfying the displayed RLLL relation for three identical Ising operators \(L^{(q)}\) is
\(R=0\).

[LEMMA] The proof is the displayed exact minor.  Any solution of all 4096 components satisfies the
32 selected equations \(B(q)r=0\).  On the localization
\(q(q-1)(q+1)\ne0\), the determinant is a unit, so \(r=0\).  Conversely, \(R=0\) satisfies the
homogeneous equation.  This proves the complete solution set in the stated graded subfamily, not
only a lower bound on its rank.

[THEOREM] For every finite ferromagnetic Ising coupling \(K>0\), one has
\(0<q=\tanh K<1\).  Therefore no nonzero—and hence no invertible—graded auxiliary \(R\) exists at
any finite ferromagnetic point in this RLLL formulation.  This conclusion was not selected or
tuned using the numerical critical-coupling benchmark.

## Gröbner bases and saturations

[COMPUTATION] All recorded characteristic-zero Gröbner calculations used graded reverse
lexicographic order over \(\mathbb Q\) with a hard 3000 s alarm.  The smallest honest fixed-point
formulation consists of the 32 selected **raw component equations**, not preassigned coordinate
equations, together with \(r_{000,000}-1\).  It completed rather than timing out.

| exact calculation at \(q=1/4\) | input formulation | result | wall time |
|---|---|---|---:|
| homogeneous ideal | 32 raw independent component forms | \(\langle r_1,\ldots,r_{32}\rangle\) | 0.169862 s |
| primary scale chart | same raw forms plus \(r_{000,000}-1\) | Gröbner basis `[1]` | 0.168703 s |
| all nonzero-coordinate charts | exactly reduced coordinate ideal plus each \(r_j-1\), 32 charts | 32/32 bases `[1]` | 0.250170 s total; 0.007419–0.008034 s each |
| invertibility saturation | reduced ideal plus \(z\det R-1\) | Gröbner basis `[1]` | 0.215650 s |

[LEMMA] The first row is an exact Gröbner reduction of the raw component ideal, so replacing those
raw generators by the 32 coordinate generators in subsequent saturation runs does not change the
ideal.  The 32 charts exhaust the complement of \(R=0\) modulo overall scalar.  Equivalently, the
projective saturation by the irrelevant coordinate ideal is the unit ideal.

[COMPUTATION] For the invertibility saturation, the generic graded \(R\) has two arbitrary
\(4\times4\) parity blocks.  Its determinant has total degree eight and 576 expanded terms.  The
Rabinowitsch equation \(z\det R-1=0\), added after the exact linear reduction, gives `[1]`.  Thus the
invertible locus is empty; this is stronger than failure of one chosen normalization coordinate.

[LEMMA] Parameter degeneracies were saturated by localizing at
\(q(q-1)(q+1)\).  The symbolic minor determinant is then a unit and forces the coordinate ideal;
combining this with any nonzero-coordinate chart or with \(z\det R-1\) gives the unit ideal.  This
localization is an exact determinant certificate rather than a separate numerical rank test.
The exceptional specializations \(q=0,+1,-1\) were deliberately removed and are not classified by
this result.

[COMPUTATION] The rational Gröbner basis completed, so neither the grade-block elimination fallback
nor the finite-field fallbacks at \(p=101\) and \(p=32003\) were run.  No modular computation is
used in the theorem.

[LEMMA] If a fallback alone had produced `[1]` for a \(\mathbb Z\)-defined saturated system modulo
one of those primes, that would prove only that the reduced ideal is the unit ideal in
\(\mathbb F_p\).  Subject to the standard good/generic-prime caveat it would be strong evidence for
characteristic-zero emptiness, but it would not replace a Gröbner basis `[1]` over \(\mathbb Q\).
The present result avoids that epistemic issue because the normalized basis over \(\mathbb Q\) is
`[1]` and the uniform-\(q\) determinant identity is also over \(\mathbb Q[q]\).

## Outcome and limitations

[THEOREM] The outcome class is **GRADED_NOGO_Q**.  The complete grade-preserving 32-coordinate
subfamily contains only the zero matrix for every \(q\notin\{0,+1,-1\}\); after scale,
nonzero-locus, or invertibility saturation its ideal is `[1]`.  There are no nonzero solutions to
classify as trivial versus interacting, and consequently no candidate whose \(2\times2\) or
\(2\times3\) transfer matrices could be tested for commutativity.

[UNRESOLVED] This theorem does not classify parity-violating \(8\times8\) matrices, auxiliary
spaces of dimension greater than two, nonidentical or independently spectral \(L\) factors,
IRF/dynamical tetrahedron relations, or singular dimension-changing projections.  It also does not
classify the exceptional algebraic points \(q=0,+1,-1\).  In particular it must not be promoted to
a no-go for every interacting auxiliary construction; it closes exactly the graded subfamily
defined above.
