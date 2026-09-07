# H-C3-34-CHAR2-EXPLICIT — fixed nonempty residual circuit construction

Owner: Main. Prospective registration; no factorization, field construction,
root extraction or numerical target verification has run. This leaves the
exhausted GF(17) factorized pool entirely; no frozen predecessor is modified.

## Fixed claim and inputs

Work in characteristic two. Fix, without search or candidate substitution,

- P(T) = T^12 + T^3 + T^2 + 1;
- f(x,y) = x^2 + y^3;
- g(x,y) = x^2 y^3 + x + y + 1;
- B = {(t^3,t^2): P(t)=0} over a splitting field K of P;
- evaluation rows x^i y^j for 0<=i<=2, 0<=j<=3.

Claim: B consists of twelve reduced affine K-rational points with both
projections injective, is exactly the proper bidegree-(2,3) complete
intersection of f and g, and has a two-dimensional vanishing pencil. Every
one-point deletion is an eleven-point circuit of this evaluation matrix;
every ten-point subset is independent. For deleted parameter a, use the
explicit remaining relation coefficients lambda_a(t)=(t+a)^2/t^4.

This is an explicit characteristic-two existence/minimality result for these
supports. It does not claim characteristic-zero or odd-characteristic
existence, classification, smallest field, global spark eleven, or absence
of smaller circuits elsewhere in the product grid. The finite GF(17)
E-empty result remains unchanged.

## Analytic proof obligation, reviewed independently before compute

1. P'=T^2 and P(0)=1, hence P has twelve distinct nonzero roots.
2. The affine intersection has no point with x*y=0; t=x/y gives the inverse
   parametrization x=t^3,y=t^2. Homogenizing f and g shows that neither x nor y
   can be infinite at a common zero. Along B the affine Jacobian determinant
   is y^2, nonzero.
3. Equal x coordinates force equal y through y=x^4+x+1; equal y coordinates
   force equal parameters because squaring is injective in a field of
   characteristic two.
4. The evaluation space becomes the ten-dimensional span of
   1,T^2,T^3,...,T^10 on the root set. The only missing degree at most ten is
   T; the duplicate T^6 and the relation T^12=T^3+T^2+1 give the two pencil
   relations f,g.
5. For a ten-root subset omitting a,b, the monic polynomial
   Q=P/((T+a)(T+b)) has coefficient of T equal to
   (a+b)/(a^2*b^2), nonzero. It cannot belong to the evaluation span, which
   proves every ten-column minor is nonsingular. Eleven columns have rank
   ten and form a circuit.
6. Lagrange interpolation on all twelve roots gives dual coefficient
   vectors t^-2 and t^-4 from the degree-eleven and degree-one coefficients.
   Their combination t^-2+a^2*t^-4 vanishes exactly at a and yields the fixed
   full-support deletion relation.

A proof-review BLOCKER stops computation until an explicit prereg/instrument
correction is committed and reviewed. A failed or unresolved proof cannot
be promoted by finite checks.

## Exact finite-field experiment

Use pinned python-flint 0.9.0 with one FLINT thread. Factor the fixed P over
F2; multiply the factors back exactly and require multiplicity one. Let m be
the least common multiple of the irreducible factor degrees, and construct
K=F_(2^m) with the first modulus returned by the library. Persist that
modulus, factorization and all element coordinates. There is no enumeration
of the 2^m field elements, no candidate search, no seed/modulus retry to
improve the outcome, and no alternative P. Factor P in K and require exactly
twelve distinct linear roots. Sort roots by their binary coefficient-vector
encoding for stable evidence order.

Run controls before this target work: a known squarefree/nonsquarefree
binary-polynomial pair, and an F4 Vandermonde with full-rank and duplicate
column controls plus a nonzero circuit relation and a corrupted relation.
They must exercise the same rank and relation procedures used for B.

Then verify and persist:

- all twelve roots, nonzero coordinates, injective x/y projections;
- f=g=0 and nonzero Jacobian determinant at all twelve points;
- full 12x12 evaluation rank exactly ten;
- all 66 distinct ten-column subsets have rank ten;
- all twelve eleven-column subsets have rank ten, and each explicit
  lambda_a vector has eleven nonzero entries and annihilates every row;
- a deliberately corrupted target relation is rejected;
- all input, source, field, matrix and result hashes.

An independent replay must import neither this instrument nor a campaign
arithmetic driver. It will use the persisted field modulus/coordinates,
independent binary-polynomial arithmetic, and recheck the point/matrix/minor
and relation claims. Source/proof review alone is not that replay.

## Resource and lifecycle contract

Each local execution uses the reviewed resource guard, one thread and nice
at least ten: sampled host CPU pause at 40%, resume below 38% for two samples,
at most two guarded jobs, at least 100 GiB free disk, sampled 2 GiB group RSS
stop, 64 MiB per-file limit, at most 8 MiB output growth and 180 seconds wall.
No priority workaround, broad search, field-element enumeration or expensive
fallback. A guard or arithmetic stop preserves all partial evidence.

Order: commit this preregistration; producer init; write and independently
review the exact proof/instrument; bind and commit final bytes; controls;
exact experiment; independent replay; owner reconciliation; one producer
freeze/close verdict; append root RESULTS and update target/progress/state.

FROZEN-CERTIFIED requires the complete analytic proof with NO BLOCKER plus
all declared exact controls, finite checks and independent replay passing.
A contradicted mathematical claim is FROZEN-NEGATIVE with a reproducer.
An instrument, resource or unresolved-proof stop is FROZEN-INCONCLUSIVE,
with the exact first failing gate and no fabricated or partial positive.
Failed attempts are retained, not overwritten or deleted.
