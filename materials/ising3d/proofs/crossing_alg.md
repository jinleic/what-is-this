# Exact local algebraic isolation of the four open-`2x3` forced crossings

Artifacts: `experiments/e96_crossing_alg.py`,
`results/spectral/crossing_alg.json`, and
`tests/test_crossing_alg.py`.

## Status and scope

**[THEOREM]** Each of the four exact c-jump brackets inherited from
`proofs/spectral_interval2.md` contains **exactly one** parameter at which

\[
f(t)=\frac{\mu_1(t)\mu_2(t)}{\mu_0(t)}
\]

is an eigenvalue of the open-`2x3` layer.  The four parameters are isolated
exact algebraic numbers: precisely, they are the `t` coordinates of unique
nonsingular real zeros of explicit systems of five polynomials over
\(\mathbb Z\), together with explicit rational isolating boxes in
`results/spectral/crossing_alg.json`.

**[THEOREM]** Those four crossings are simple (the corresponding five
hypersurfaces meet transversely), and there is no second forced crossing in
any of the four inherited fine brackets.

**[COMPUTATION]** The exact sector characteristic polynomials, their
factorizations, the Krawczyk interval inequalities, the inertia counts, and
the derivative bounds are finite calculations rebuilt by the standalone test.

**[UNRESOLVED]** This is *bracket-complete*, not scan-complete.  No
minimal polynomial in \(\mathbb Z[t]\), global elimination polynomial
\(R(t)\), or Sturm count excluding further crossings in the unbracketed parts
of \([1/4,7/20]\) is claimed.  The direct first nested resultant expansion
recorded below reached its 1,200-second wall without producing a polynomial.

The exact brackets and the ordered fourth branch are:

| \(r\) | exact bracket \(B_r\) | branch equal to \(f\) | sector of that branch | Krawczyk norm upper bound |
|---:|---|---|---:|---:|
| 1 | \([2248027753/8388608000,4496055507/16777216000]\) | \(\mu_4\) | 4 | \(1/1000\) |
| 2 | \([4881706207/16777216000,152553319/524288000]\) | \(\mu_5\) | 5 | \(1/1000\) |
| 3 | \([2522992499/8388608000,5045984999/16777216000]\) | \(\mu_6\) | 6 | \(1/1000\) |
| 4 | \([553214907/1677721600,5532149071/16777216000]\) | \(\mu_7\) | 4 | \(1/500\) |

The decimal guides in the JSON are orientation only and have no role in any
final inequality.

## 1. Exact polynomial reduction

Let \(M(t)=(2t)^7q(t)^4R(t)\) be the positive integer-polynomial
representative from the wave-9 note, with
\(q(t)=(1+t^2)/(2t)\).  Direct expansion over all 64 layer states gives a
symmetric matrix with nonnegative integer coefficients and degree at most 26.
The same direct expansion proves the reciprocal identity

\[
M(t)=t^{26}M(t^{-1}).
\]

**[LEMMA — reciprocal sector reduction.]** For every joint character sector
of row reflection, column reflection, and global spin flip,
\(t^{-13}M(t)\) is a matrix over \(\mathbb Z[z]\), where
\(z=t+t^{-1}\).

*Proof.*  Each scalar entry is a palindromic degree-26 polynomial.  Writing
its coefficient pair at powers \(13\pm k\) as a multiple of
\(t^k+t^{-k}\), and using
\(C_0=2,C_1=z,C_{k+1}=zC_k-C_{k-1}\), gives a polynomial in \(z\).
Integer character vectors preserve that property. \(\square\)

The independently reconstructed sector dimensions are

\[
(14,6,6,6,10,10,6,6).
\]

For a sector \(s\), let \(G_s\) be its diagonal integer Gram matrix and
\(B_s(z)\) its exact character-basis block.  The producer factors

\[
\chi_s(X,z)=\det(XG_s-B_s(z))\in\mathbb Z[X,z].
\]

The low-branch factors \(Q_s\) needed here, and their primitive restored
polynomials,

\[
P_s(x,t)=\operatorname{prim}\left(t^{13\deg_XQ_s}
  Q_s(x/t^{13},t+t^{-1})\right)\in\mathbb Z[x,t],
\]

have the following exact bidegrees.  The SHA-256 digests in the JSON bind the
coefficient lists that the clean-room verifier reconstructs.

| sector \(s\) | \(\deg_XQ_s\) | \(\deg_zQ_s\) | \(\deg_tP_s\) | terms of \(P_s\) |
|---:|---:|---:|---:|---:|
| 0 | 14 | 133 | 315 | 2040 |
| 4 | 9 | 90 | 207 | 880 |
| 5 | 9 | 81 | 198 | 797 |
| 6 | 6 | 59 | 137 | 382 |
| 7 | 6 | 55 | 133 | 356 |

Sectors 4 and 5 each have one additional explicit linear factor; exact
inertia puts its root outside every low-eigenvalue isolation window used
below.  Thus the selected low roots are roots of the displayed \(P_s\), not
of an omitted linear factor.

## 2. Five-polynomial certificates

For crossing \(r\), set \(d_r=4,5,6,7\) and
\(s_r=4,5,6,4\), respectively.  In variables \((t,a,b,c,d)\), define

\[
\mathcal F_r=
\big(P_0(a,t),\;P_5(b,t),\;P_7(c,t),\;P_{s_r}(d,t),\;bc-ad\big).
\]

Every coefficient of every member of \(\mathcal F_r\) is an exact integer.
The last equation is exactly the scale-invariant forced-value relation
\(d=bc/a\) because \(a>0\).

For each \(r\), the JSON stores a rational center \(c_r\), radii

\[
\left(5\cdot10^{-11},\;2\cdot10^{-10},\;2\cdot10^{-10},
2\cdot10^{-10},\;2\cdot10^{-10}\right),
\]

and a rational 5-by-5 preconditioner \(C_r\).  Let \(D_r\) be that closed
box.  The verifier evaluates the polynomial Jacobian at the exact center,
bounds its variation on \(D_r\) by coefficientwise absolute Hessians, and
checks the interval Krawczyk inclusion

\[
K_r(D_r)\subset\operatorname{int}D_r,
\qquad
\|K_r-c_r\|_{D_r,\infty}<
\begin{cases}
1/1000,&r=1,2,3,\\
1/500,&r=4.
\end{cases}
\]

(The coordinatewise inclusion ratio is also below \(1/100\) in every case.)
All operations in this test are `Fraction` operations after the stored decimal
preconditioner has been read as an exact rational number.

By the Krawczyk theorem, each \(\mathcal F_r\) has exactly one real zero in
\(D_r\), and the nonzero-Jacobian condition holds there.  This gives a
zero-dimensional rational polynomial system locally, so every coordinate of
that zero, in particular its \(t\) coordinate, is algebraic.  It also proves
transversality: eliminating the four simple sector-root directions from the
Jacobian leaves the derivative of \(bc-ad\) along the four branches, which
cannot vanish when the five-by-five Jacobian is nonsingular.

## 3. Identification with the ordered spectrum

The Krawczyk calculation by itself isolates a product coincidence among
sector roots.  The following exact checks identify it with the forced value
of the *ordered* spectrum.

At each rational `t` center, fraction-free Sylvester inertia of the eight
character-sector pencils gives the threshold count vectors stored in the JSON.
The five thresholds separate, in order,

\[
\mu_0,\quad\mu_1,\quad\mu_2,\quad\mu_{d_r}.
\]

At the same center, an exact \(10^{-10}\) window around each of the four
selected roots contains exactly one eigenvalue, in the stated sector, and no
other sector eigenvalue.  The test repeats every one of those integer-pivot
inertia calculations.

For transport from the center to the inherited fine bracket \(B_r\), exact
coefficient bounds give

\[
\|M'(t)\|_2\,5\cdot10^{-11}<10^{-7}.
\]

The rational `t` center itself lies in \(B_r\), and all terms used in the
norm bound have nonnegative coefficients, so evaluation at the rational upper
endpoint of \(B_r\) covers the whole center-to-bracket path.  The analogous
selected-sector bounds hold as well.  Exact inertia at every separator plus
and minus \(10^{-7}\) shows that none of these spectral separators can be
crossed anywhere in \(B_r\).  Therefore the labels

\[
a=\mu_0,\qquad b=\mu_1,\qquad c=\mu_2,\qquad d=\mu_{d_r}
\]

are stable throughout the inherited bracket.

Finally, exact implicit differentiation of each \(P_s(x,t)=0\) on the
larger \(5\cdot10^{-9}\) eigenvalue tube yields

\[
|\lambda'(t)|<1/20.
\]

The lower bound for \(|P_{s,x}|\) used in that quotient is positive on the
entire tube.  A first-exit argument therefore bootstraps the local implicit
branch across \(B_r\): a branch cannot leave the tube before the displayed
derivative bound rules out its leaving it.  Starting from the exact
\(10^{-10}\) center enclosure, every selected branch remains within

\[
\delta=10^{-10}+(1/20)(5\cdot10^{-11})
 =\frac{41}{400000000000}<2\cdot10^{-10}
\]

of its stored center throughout \(B_r\).

This last estimate supplies the missing converse from a generic forced
equality to the selected branch.  Let \(a_c,b_c,c_c\) denote the exact
centers.  Positivity and the branch bounds give, at every \(t\in B_r\),

\[
\frac{(b_c-\delta)(c_c-\delta)}{a_c+\delta}
\ \leq\ f(t)\ \leq\
\frac{(b_c+\delta)(c_c+\delta)}{a_c-\delta}.
\]

The two rational endpoints are computed and stored for every bracket.  Exact
comparison gives the strict inclusion

\[
q_{d,-}<\frac{(b_c-\delta)(c_c-\delta)}{a_c+\delta}
\leq\frac{(b_c+\delta)(c_c+\delta)}{a_c-\delta}<q_{d,+},
\]

where \(q_{d,-},q_{d,+}\) are the final two stored separators.  Their
protected inertia counts differ by exactly one.  Thus **every** forced
equality in \(B_r\), not merely an equality suggested by the c-jump, is with
the sole eigenvalue in \((q_{d,-},q_{d,+})\), namely \(\mu_{d_r}=d\).
All four participating ordered branches lie in their Krawczyk eigenvalue
boxes, so every forced crossing in \(B_r\) yields a zero of
\(\mathcal F_r\) in \(D_r\).

The prior theorem in `proofs/spectral_interval2.md`, independently backed by
its c-jump endpoint certificates, provides at least one forced crossing in
each \(B_r\), with exact c-count jumps
\(4\to5,5\to6,6\to7,7\to8\).  Each \(B_r\) lies in the corresponding
Krawczyk `t` box.  Krawczyk uniqueness therefore proves both existence and
uniqueness in that original rational bracket.  In particular, the unique
Krawczyk zero is the inherited forced crossing rather than a merely formal
sector-product solution.

## 4. Resultant attempt and honest remainder

A direct global approach would eliminate the four sector eigenvalue variables
from the preceding equations, then factor the resulting univariate
\(R(t)\) and use exact Sturm counts on \([1/4,7/20]\).  The reciprocal form
was used to reduce the first stage to

\[
U(z,y)=\operatorname{Res}_b\left(Q_5(b,z),\;b^6Q_7(y/b,z)\right),
\]

whose inputs have \(X\)-degrees \((9,6)\) and term counts \((309,134)\).

**[COMPUTATION — wall certificate.]** An exact SymPy multivariate resultant
attempt was allowed 1,200 seconds.  It did not produce \(U(z,y)\); later
eliminations and the putative global \(R\) were consequently not attempted.
This wall is recorded as data, not promoted into a negative theorem.

Accordingly, the following remain **[UNRESOLVED]**:

1. minimal polynomials and conventional one-variable isolating intervals for
   the four \(t_r^*\);
2. a global resultant/Sturm proof that the four brackets exhaust all forced
   crossings in \([1/4,7/20]\);
3. the presence or absence of any additional crossings in the five
   unbracketed scan runs.

## 5. Standalone verification

Run, from the repository root:

```bash
timeout 1800 .venv/bin/python tests/test_crossing_alg.py
```

The test imports neither the producer nor a generated polynomial.  It rebuilds
from raw 2-by-3 states and exact integer arithmetic, checks the low-sector
factor digests, verifies every rational Krawczyk inclusion, recomputes the
center/margin inertia counts and branch derivative bounds, and checks that the
four brackets and c-count jumps agree with the independently tested wave-9
artifact.
