# Exact finite certificates for three $K_c$ barrier routes

Artifacts: `experiments/e85_kc_interval2.py`,
`results/bounds/kc_interval2.json`, and `tests/test_kc_interval2.py`.

## Result

**[UNRESOLVED]** No route in this audit changes the certified simple-cubic
interval:

\[
  0.2122119011661678393310862783954278184914
  \le K_c \le
  0.2527310098586630030260020266135701299926.
\]

The lower endpoint remains
\(\operatorname{atanh}(c_{36}^{-1/36})\), with

\[
 c_{36}=2941370856334701726560670.
\]

This note adds three independently checkable pieces of finite/exact progress:

1. **[COMPUTATION]** a fixed family of four Simon--Lieb prisms with 18- and
   20-site cross-sections is evaluated by the full exact parity transfer;
   all four reject the exact challenge point above the incumbent lower bound.
2. **[LEMMA]** an exact finite-torus infrared/GKS linear-program chain gives
   positive finite-volume zero-mode floors at both stored sides $L=4,6$;
   notably it gives a positive $L=6$ floor at $K=1/5$ where the naïve
   all-ceiling bound is negative.
3. **[THEOREM]** a first-backtrack injection gives a strictly sharpened,
   FKG-free self-avoiding-walk submultiplicativity inequality for every
   $m,n\ge1$.  Its available exact numerical consequences do not beat the
   incumbent, and the exact missing $n=30$ ratio lemma is isolated.

All comparison decisions in the artifact are integer or `Fraction`
comparisons.  The 90-dps interval evaluation of `atanh` is formatting only.
No critical-coupling benchmark is used to select a shape, point, or result.

## 1. Wider exact Simon--Lieb transfer certificates

For a free rectangular box $B$, let $P_B(t)$ be its even-subgraph polynomial
and let $Q_B(t)$ be the exact boundary-weighted polynomial used in the audited
Simon--Lieb criterion.  The criterion is

\[
 \kappa_B(t)=\frac{tQ_B(t)}{P_B(t)}<1
 \quad\Longrightarrow\quad K_c\ge\operatorname{atanh}(t).
\]

**[COMPUTATION]** The exact transfer evaluates the sign of
\(tQ_B(t)-P_B(t)\) as an integer after clearing the rational denominator.  It
uses all $2^{ab}$ transfer states for an $a\times b$ cross-section.  The
certificate family was fixed as

\[
 (3,6,4),\quad(3,6,12),\quad(4,5,4),\quad(4,5,16),
\]

before its non-authoritative binary64 locators were run.  Thus every member has
cross-section $ab>16$, strictly beyond the previous audited family.

| box $B$ | $ab$ | exact safe $t$ | outward lower endpoint for $\operatorname{atanh}(t)$ | sign at safe $t$ | sign at $23/110$ |
|---|---:|---:|---:|---:|---:|
| $3\times6\times4$ | 18 | $1993/10000$ | $0.2020034935795436133893256675929987185959$ | $-$ | $+$ |
| $3\times6\times12$ | 18 | $2023/10000$ | $0.2051295405318898002935737421801965685413$ | $-$ | $+$ |
| $4\times5\times4$ | 20 | $1003/5000$ | $0.2033576322702597894725981993171751197084$ | $-$ | $+$ |
| $4\times5\times16$ | 20 | $2039/10000$ | $0.2067983799232241126594147841719378985339$ | $-$ | $+$ |

Here a minus sign means the exact integer $tQ_B-P_B$ is negative.  The JSON
records each signed integer's bit length and SHA-256 of its unsigned
big-endian encoding, so the displayed compact table remains independently
addressable without placing multi-thousand-bit integers in this note.

**[COMPUTATION]** The strongest safe point is the $4\times5\times16$ point.
The two exact power comparisons

\[
  2039^{36}c_{36}<10000^{36},\qquad
  23^{36}c_{36}>110^{36}
\]

show respectively that this entire fixed family remains below the incumbent
lower endpoint and that the common challenge $23/110$ lies strictly above it.
The four positive challenge residuals rule out an improvement by this completed
family.  This is a finite-family obstruction only: it makes no claim about
larger cross-sections, other geometries, or sharper Simon--Lieb contractions.

The producer also has two exact controls: complete $2\times2\times2$
spin enumeration agrees with the parity transfer, and on $2\times3\times2$
the streamed residual equals an independently reconstructed $pQ-qP$.

## 2. Exact finite-volume infrared/GKS floor

Let $T_L=(\mathbb Z/L\mathbb Z)^3$, $N=L^3$, and let $x_j=\widehat G(k)$ be
one common Fourier value on a nonzero signed-permutation orbit $O_j$.  Write
$m_j=|O_j|$, $\lambda_j=3-\sum_i\cos k_i$, and

\[
 A_{zj}=\sum_{k\in O_j}\cos(k\mathbin\cdot z),\qquad
 S=\frac1N\sum_jm_jx_j,\qquad M_L^2=1-S.
\]

$M_L^2$ is the finite-volume zero-mode mass used by the infrared argument; it
is not an assertion of a nonzero one-point magnetization on a finite
zero-field torus.

**[LEMMA]** Under exactly the finite-torus constraint system audited in
`proofs/upper_infrared.md`, the following rational linear program is valid:

\[
\begin{aligned}
0&\le x_j\le\frac1{2K\lambda_j},\\
\sum_j(m_j-A_{zj})x_j&\le N &&\text{for every }z\ne0,\\
\sum_j(A_{zj}-m_j)x_j&\le0 &&\text{for every }z\ne0,\\
\sum_jm_jx_j&\le N,\\
\sum_j(m_j-A_{ej})x_j&\le\frac{N-1}{6K}.
\end{aligned}
\tag{1}
\]

Indeed, the first line is the infrared ceiling.  The next two are respectively
$G_L(z)\ge0$ and $G_L(z)\le1$ after using $G_L(0)=1$; the fourth is
$M_L^2\ge0$; and the last is the finite-volume nearest-neighbour energy
constraint.  Maximizing $\sum_jm_jx_j$ in (1) therefore yields the rigorous
finite-volume conclusion

\[
 M_L^2\ge1-S_*(K,L).\tag{2}
\]

**[COMPUTATION]** Exact `Fraction` primal and dual solutions have equal
objectives for every entry below; the full vectors, constraint labels, and
bases are stored in the JSON and rebuilt by the standalone verifier.

| $K$ | $1-S_*(K,4)$ | $1-S_*(K,6)$ |
|---:|---:|---:|
| $1/5$ | $251/4608$ | $253/15552$ |
| $21/100$ | $607/8064$ | $133/5832$ |
| $11/50$ | $863/8448$ | $34663/940896$ |
| $6/25$ | $1631/9216$ | $737509/7185024$ |

The strengthening at $L=6,K=1/5$ is concrete.  The stored exact Green data are

\[
 C_6(0)=\frac{1289503}{2993760},
 \qquad \min_z C_6(z)=-\frac{13861}{598752},
\]

whereas the naïve all-ceiling value gives only

\[
 1-\frac{C_6(0)}{2(1/5)}=-\frac{91999}{1197504}<0.
\]

Equation (2), using the exact LP optimum, instead gives
\(M_6^2\ge253/15552>0\).  This is an exact finite-volume improvement, not a
thermodynamic conclusion.

For a useful diagnostic, the literal all-ceiling profile becomes feasible only
at or above

\[
 \frac{C_4(0)-\min C_4}{2}=\frac{13}{60},\qquad
 \frac{C_6(0)-\min C_6}{2}=\frac{5147}{22680},
\]

and is directly checked feasible at $K=6/25$ in both systems.  Thus the
finite LP is genuinely tighter before that profile is feasible.

**[UNRESOLVED] Missing uniform lemma.** To turn (2) into a lower bound on the
thermodynamic $K_c$, it would suffice to prove that there exist
$K<I_3/2$, $m>0$, and $L_0$ such that

\[
 M_L^2(K)\ge m\quad\text{for every even }L\ge L_0.\tag{3}
\]

The displayed certificates establish (2) separately for $L=4,6$ and cannot
supply the quantifiers in (3).  The finite checks therefore do not evade the
uniform-in-$L$ barrier identified in `proofs/upper_infrared.md`.

## 3. FKG-free self-avoiding-walk refinement

Let $c_n$ be the number of $n$-step cubic self-avoiding walks from the origin,
and let $a_r$ count the $r$-step walks which avoid one specified neighbour of
the origin.

**[THEOREM] First-backtrack injection.** For every $m,n\ge1$,

\[
 c_{m+n}\le c_m\bigl(c_n-a_{n-1}\bigr).\tag{4}
\]

*Proof.*  Choose an $m$-step walk $\omega$.  For every $(n-1)$-step walk
$\rho$ from its penultimate vertex which avoids its endpoint, prepend to
$\rho$ the reverse of $\omega$'s final edge.  The resulting $n$-step
walk $\eta$ is self-avoiding, but its concatenation with $\omega$ revisits
$\omega_{m-1}$ and is invalid.  Translation gives exactly $a_{n-1}$ choices
per $\omega$, and the construction is injective: $(\omega,\eta)$ recovers
$\rho$ after removing the forced first edge.  Since valid concatenated pairs
are in bijection with $(m+n)$-step walks, removing these invalid pairs proves
(4). $\square$

**[COMPUTATION]** Direct visited-set backtracking gives

\[
(c_0,\ldots,c_9)=(1,6,30,150,726,3534,16926,81390,387966,1853886),
\]
\[
(a_0,\ldots,a_8)=(1,5,25,121,589,2821,13565,64661,308981).
\]

These data check (4) for all $m+n\le9$.  With the audited external counts
$c_{30},\ldots,c_{36}$, they also check it in the six mixed cases
$(m,n)=(30,6),\ldots,(35,1)$.

**[EXTERNAL]** The source for $c_{30},\ldots,c_{36}$ is the cached
`schram_barkema_bisseling2011.pdf`, whose SHA-256 is

```text
898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12
```

The test hashes that exact file before accepting any statement which uses those
counts.

Fekete submultiplicativity gives
\(\mu\le c_{36}^{1/36}\).  Consequently the exact integer inequality

\[
 c_{36}c_n^{18}\le c_{n+2}^{18}\tag{5}
\]

implies \(\mu^2\le c_{n+2}/c_n\).  **[COMPUTATION]** (5) holds for
$n=0,1,\ldots,6$ and fails at $n=7$; all signs and hashes are in the JSON.
This is a valid finite-ratio consequence of Fekete plus the external exact
count, but it is not the proposed general Kesten finite-ratio theorem.

**[COMPUTATION]** Applying (4) to every available external pair
$30\le m\le36$, $1\le n\le8$ yields no self-avoiding-walk growth bound
stronger than $c_{36}^{1/36}$.  Thus the injection sharpens finite
submultiplicativity but does not change the lower endpoint.

**[UNRESOLVED] Exact improving lemma.** The exact comparison

\[
 c_{32}^{18}<c_{30}^{18}c_{36}
\]

shows that a first-principles proof of
\(\mu^2\le c_{32}/c_{30}\) would imply
\(v_c\ge\sqrt{c_{30}/c_{32}}>c_{36}^{-1/36}\), strictly improving the
current lower endpoint.  No such $n=30$ finite-ratio proof is supplied here.

## Scope and independent verification

**[COMPUTATION]** Run

```sh
timeout 3600 .venv/bin/python experiments/e85_kc_interval2.py
timeout 3600 .venv/bin/python tests/test_kc_interval2.py
```

from the repository root.  The producer writes the exact certificate envelope;
the standalone verifier never imports the producer.  It independently rebuilds
full-transfer route-1 residuals, all $L=4,6$ rational LP constraints and
primal-dual equalities, the small self-avoiding-walk data, every displayed
finite integer comparison, and the external-source SHA-256.

**[COMPUTATION]** The recorded producer elapsed time was
`1332.731340790997` seconds; the independently timed standalone verifier took
`1018.04` seconds.  Both commands used a 3600-second timeout, exceeding twice
the observed wall time.
