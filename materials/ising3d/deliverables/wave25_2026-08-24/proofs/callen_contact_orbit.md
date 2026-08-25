# Contact-inclusive Callen orbit obstruction on the periodic cubic cell

Artifacts: `experiments/e250_callen_contact_orbit.py`,
`results/correlations/callen_contact_orbit.json`, and the clean-room verifier
`tests/test_callen_contact_orbit.py`.

## 1. Theorem and exact scope

**[THEOREM].** Let `0<v=tanh(K)<1` and consider the zero-field nearest-neighbour Ising model on the
simple six-regular periodic graph `C3 square C3 square C3`. Take every Callen identity whose pivot is
any vertex and whose pivot-independent factor is `sigma_B` for an odd subset `B` of the six pivot
neighbours. Quotient these `27*32=864` rows and their correlation columns by translations and the
full signed-coordinate-permutation cubic group. No nonzero linear combination of the four resulting
row classes can cancel every column except normalization and the nearest-neighbour pair orbit.

The tested block includes marked-spin contacts, all pivots, both three-neighbour geometries, and the
size-five neighbour-contact reduction of the `U=N` row excluded from e246's far-mark matrix. It
therefore crosses the named pair-row boundary of e246 rather than repeating it.

**[SCOPE].** This is a linear, radius-one, cubic-symmetry-aggregated obstruction. It does not cover
radius-expanded factors, nonlinear identities, auxiliary observables, or nonsymmetric row
combinations with zero cubic average. It is not a susceptibility closure theorem and does not solve
the finite or infinite three-dimensional model.

## 2. The 864 contact rows reduce to four classes

Write the pivot as `0`, its six neighbours as the signed coordinate directions, and

\[
 \tanh\!\left(K\sum_{i\sim0}\sigma_i\right)
 =c_1(v)e_1+c_3(v)e_3+c_5(v)e_5,
\]

where `e_j` is the elementary symmetric polynomial in the neighbour spins. For each odd
`B subset N(0)`, the exact row is

\[
 C(\{0\}\cup B)-
 \sum_{\substack{A\subset N(0)\\ |A|\text{ odd}}}
 c_{|A|}(v)C(B\mathbin\triangle A)=0.                 \tag{2.1}
\]

The origin stabilizer has order `48`: coordinate permutations and independent sign reversals.
Its odd-subset orbits have direction-count types

\[
 (1,0,0),\qquad (2,1,0),\qquad (1,1,1),\qquad (2,2,1),
\]

of sizes `6,12,8,6`. Translation over 27 pivots therefore gives raw row-orbit sizes

\[
 162,\quad324,\quad216,\quad162,
\]

which sum to all `864` rows. The producer derives these orbits directly from all 1,296 space-group
permutations; the verifier independently constructs `S3 wr S3`.

Canonicalizing every support in (2.1) yields nine correlation-column orbits: normalization, the
nearest-neighbour pair, one other pair orbit, four four-spin orbits, and two six-spin orbits. Remove
the first two allowed columns. The remaining block is `4 x 7`.

## 3. A four-column minor is positive at every physical coupling

Exact elimination at the rational control `v=1/3` selects four nonallowed columns. Keeping the
formal variables `c1,c3,c5`, their determinant factors without approximation as

\[
 \Delta_c=4(c_1-c_5)(c_1+4c_3+c_5).                    \tag{3.1}
\]

The exact Callen coefficients inherited from e240 are

\[
\begin{aligned}
 c_1&=\frac{v(v^8+16v^6+46v^4+16v^2+1)}{D(v)},\\
 c_3&=\frac{-2v^3}{(1+v^2)(1+14v^2+v^4)},\\
 c_5&=\frac{16v^5}{D(v)},\\
 D(v)&=(1+v^2)(1+6v^2+v^4)(1+14v^2+v^4).
\end{aligned}                                           \tag{3.2}
\]

Substitution gives

\[
\boxed{
 \Delta(v)=
 \frac{4v^2
 (v^8+16v^6+30v^4+16v^2+1)
 (v^8+8v^6+14v^4+8v^2+1)}
 {(1+v^2)^2(1+6v^2+v^4)^2(1+14v^2+v^4)^2}.}             \tag{3.3}
\]

Every displayed numerator and denominator factor is strictly positive for `v>0`. Hence
`Delta(v)>0` throughout `0<v<1`, so the nonallowed `4 x 7` block has full row rank four. Its left
kernel is zero. A row combination cancelling all seven nonallowed correlation orbits must therefore
be the zero combination and cannot retain a nearest-neighbour pair coefficient. This proves the
theorem.

## 4. Verification boundary

The producer derives the complete space-group and support-orbit tables, checks all 864 raw rows are
covered, reconstructs the e240 coefficient functions, verifies the `v=1/3` rank, factors (3.1) and
(3.3), and fails before writing on any semantic, provenance, CPU, or RSS gate. The verifier imports
neither e250 nor e240: it rebuilds the group as `S3 wr S3`, reads only the stored e240 coefficient
record, reconstructs the literal four-row minor, and rechecks every factor and source hash.
