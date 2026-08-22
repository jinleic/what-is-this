# Exact specialization locus for the one-line-cabled \(16\times16\) tetrahedron system

[DEFINITION] This note uses exactly the fixed-order one-line-cabled auxiliary-\(R\) system of `proofs/tetra16.md`: its coefficient matrix is \(B(q)\in\mathbb Z[q]^{16384\times256}\), split into two 128-column parity sectors.  It does not concern a free \(\mathbb C^4\)-leg local operator or solve the three-dimensional Ising model.

## Result

[THEOREM] Over every characteristic-zero field containing the specialization \(q\),

\[
 \operatorname{rank}B(q)<256\quad\Longleftrightarrow\quad q\in\{-1,0,1\}.
\]

The exact exceptional ranks are

| \(q\) | preserving-sector rank | reversing-sector rank | \(\operatorname{rank}B(q)\) | nullity |
|---:|---:|---:|---:|---:|
| \(-1\) | 0 | 0 | 0 | 256 |
| \(0\) | 14 | 16 | 30 | 226 |
| \(1\) | 29 | 29 | 58 | 198 |

[COMPUTATION] Each displayed rational specialization rank is the exact SymPy \(\mathbb Q\)-rank of every one of the 8192 rows in its sector, not a modular estimate.

## Two complementary full minors

[DEFINITION] Let \(P_0,P_1\in\mathbb Z[q]\) be the preserving and reversing 128-by-128 minors using the inherited primary row sets, and let \(A_0,A_1\) use the disjoint reverse-scan row sets.  Define

\[
 D_P=P_0P_1,\qquad D_A=A_0A_1.
\]

[COMPUTATION] All four polynomials were reconstructed without forming a dense symbolic determinant.  For each minor, its entry valuations and degrees give a support-matching bound; its shifted determinant was evaluated at every integer node \(1,\ldots,D+1\) over 17 exact primes, Newton-interpolated over each finite field, then coefficientwise CRT-lifted past a proven height bound.  The resulting data are:

| minor | valuation lower bound | degree bound | actual degree | shifted degree | coefficient hash |
|---|---:|---:|---:|---:|---|
| \(P_0\) | 283 | 871 | 868 | 585 | `a160d06753a6d62e6c29330c2cbac2aeb21d7b745cf3a58b2867ceef94c4d1b1` |
| \(P_1\) | 271 | 880 | 872 | 601 | `4d6b9a3bbbd78568f964a533e6789416d2ccb05f649e61787aff99eb83c47369` |
| \(A_0\) | 665 | 1253 | 1240 | 575 | `ef2cb662801b141d262c5d38838c88d6f0e8de29b7c8632fdaf33c625d56e26a` |
| \(A_1\) | 656 | 1265 | 1248 | 592 | `334fc6dd997f240578b434421023efc66f4f0260f78fb6a9c47a00b588f3229d` |

[COMPUTATION] The maximum shifted-degree envelope is 609.  Each CRT modulus has 526 bits, while the independently derived coefficient-height bounds have at most 512 bits.  Thus the stored coefficient vectors in `results/integrability/tetra16_locus.json` are uniquely determined integer vectors, not numerical fits.  The producer also performs exact Bareiss checks at \(q=-3,3/7\); the clean-room verifier uses the distinct points \(q=-5,5/7\).

[THEOREM] The exact complementary-minor gcd is

\[
 \gcd_{\mathbb Q[q]}(D_P,D_A)
 =G(q)=q^{584}(q-1)^{580}(q+1)^{546}(q^2+1)^4. \tag{1}
\]

[COMPUTATION] An exact \(\mathbb Z[q]\) gcd completed in 6.980353 seconds.  It has degree 1718 (coefficient-vector SHA-256 `048ff238f13fda7660865cd03ce28119d787ec4b53d0e3c8ab8f8117d2a53ac3`); the exact quotient degrees are 22 for \(D_P/G\) and 770 for \(D_A/G\).  The factorization of \(G\) completed in 21.511872 seconds.

[THEOREM] Equation (1) is independently certified by the explicit residue-field Bézout identity stored in the result artifact.  After division by \(G\), the two residual polynomials have gcd 1 over \(\mathbb F_p[q]\) for

\[
 p=1{,}999{,}999{,}973,
 \qquad \operatorname{LC}(D_P/G)\equiv1{,}999{,}950{,}821\pmod p.
\]

The artifact contains the exact low-to-high coefficient lists \(u,v\) and the standalone verifier checks

\[
 u(D_P/G)+v(D_A/G)=1\pmod p.
\]

If a nonconstant primitive rational common divisor remained, its leading coefficient would divide that of \(D_P/G\), hence would remain nonzero in this selected residue field and yield a nonconstant modular gcd—a contradiction.  This provides the exact upper bound matching the direct gcd computation.

## From the minor gcd to the full system

[LEMMA] If \(\operatorname{rank}B(q)<256\), then \(G(q)=0\).

[THEOREM] Both \(D_P\) and \(D_A\) are 256-by-256 block-diagonal minors of \(B(q)\).  A rank drop forces both to vanish, hence forces their gcd to vanish.  Therefore (1) proves \(\operatorname{rank}B(q)=256\) for every \(q\notin\{0,1,-1,i,-i\}\).

[THEOREM] The two remaining roots of (1) are not exceptional.  For each \(q=\pm i\), two split-prime reductions contain nonzero recorded 128-by-128 minors in **both** parity sectors:

| algebraic point | reduction \(p,q\bmod p\) | preserving determinant | reversing determinant |
|---|---|---:|---:|
| \(i\) | \(5,2\) | 4 | 4 |
| \(i\) | \(13,5\) | 10 | 10 |
| \(-i\) | \(5,3\) | 1 | 1 |
| \(-i\) | \(13,8\) | 3 | 3 |

For example, \(i\mapsto2\in\mathbb F_5\) and \(i\mapsto5\in\mathbb F_{13}\), while \(2^2+1\equiv0\pmod5\) and \(5^2+1\equiv0\pmod{13}\).  Each nonzero 128-by-128 minor gives sector rank at least 128 over \(\mathbb Q(i)\); the exact 128-column upper bound gives equality.  Hence both sectors, and therefore the full system, have rank 256 at \(q=i\) and \(q=-i\).  The two-prime certificates are recorded with their row sets in the artifact and independently rebuilt by the standalone verifier.

Combining this result with the exact rational ranks in the first table proves the stated complete locus.

## Reproduction

[CHECK] Regenerate the exact certificate (measured wall: 674.57 seconds):

```bash
timeout 1800 .venv/bin/python experiments/e95_tetra16_locus.py
```

[CHECK] Run the clean-room standalone verifier (measured wall: 747.64 seconds):

```bash
timeout 1800 .venv/bin/python tests/test_tetra16_locus.py
```

The verifier does not import the producer.  It rebuilds the seven-bit operators, all four minors, enough modular interpolation nodes and CRT bounds, fresh exact determinant values, the final Bézout identity, all three rational rank drops, and the four split-prime algebraic full-rank minors.

[SCOPE] This is an exact all-specialization theorem only for the fixed-order one-line-cabled 16-by-16 auxiliary-\(R\) ansatz.  It does not establish an RLLL solution for an unrestricted \(\mathbb C^4\) leg, cabling of all three auxiliary lines, spectral/nonidentical local factors, or the 3D Ising model.
