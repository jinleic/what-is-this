# Gate B proof: the repair ratio is unbounded

## Theorem

Let \(A_+\), \(\varepsilon_\vee\), the cap, and the Reimer incidence
condition be exactly as frozen in [DEFINITIONS.md](DEFINITIONS.md). Then

\[
\boxed{c_{\rm cl}^\star=+\infty.}
\]

More explicitly, there is a normalized cap/Reimer family
\(\mathcal B\subseteq2^{[7]}\) such that, for every integer \(k\ge1\),
\(\mathcal F_k=\mathcal B^{\boxtimes k}\subseteq2^{[7k]}\) is normalized and
cap/Reimer, and

\[
\varepsilon_\vee(\mathcal F_k)=1-\left(\frac{17}{81}\right)^k,
\qquad
A_+(\mathcal F_k)=kA_+(\mathcal B),
\]

with the certified bound

\[
A_+(\mathcal B)<-\frac7{250}.
\]

Consequently

\[
\frac{-A_+(\mathcal F_k)}{
\varepsilon_\vee(\mathcal F_k)}
>
\frac{7k/250}{1-(17/81)^k}
>
\frac{7k}{250}\longrightarrow+\infty.
\]

Thus no fixed scalar coefficient \(c\) can make
\(A_+(\mathcal F)+c\varepsilon_\vee(\mathcal F)\ge0\) on every admissible
family.

## Certified base lemma

Let \(\mathcal B\) be the 45-row family in
[`candidates/n7_block_extremizer.json`](candidates/n7_block_extremizer.json).
It is the union of the cells

\[
(0,1),(0,3),(1,0),(1,2),(1,5),(2,0),(2,1)
\]

for the coordinate partition \(\{5,6\}\sqcup\{0,1,2,3,4\}\). The standalone
checker [`verify_gate_b.py`](verify_gate_b.py), which imports no other
`math/uc` module, proves the following.

1. \(|\mathcal B|=45\), every coordinate occurs 18 times, and the total
   incidence is 126.
2. The exact Reimer threshold is 124. In particular, both cap and Reimer hold.
3. The seven columns are active and pairwise distinct.
4. Exactly 1,600 of the \(45^2\) ordered join pairs leave the family, so
   \(\varepsilon_\vee(\mathcal B)=64/81\) and the join-success probability is
   \(17/81\).
5. There are 240 coordinate automorphisms and 21 equal-size order orbits.
6. A 256-bit Arb recurrence gives a rigorous upper certificate for every
   one-sided Bellman maximum. At each of the 982 clamp-ambiguous states it
   retains the unconstrained concave maximum, which upper-bounds the constrained
   maximum, and gives

   \[
   A_+(\mathcal B)
   \le -0.0286491867944683178160269789255044875\ldots
   <-\frac7{250}.
   \]

The authoritative interval artifact is
[`certificates/gate_b_unbounded_arb_v4.json`](certificates/gate_b_unbounded_arb_v4.json);
versions 1--3 are preserved intermediate audit checkpoints.
The verifier explicitly constructs the 2,025-row square
\(\mathcal B^{\boxtimes2}\) and counts 3,920,000 missing ordered joins, giving

\[
\varepsilon_\vee(\mathcal B^{\boxtimes2})
=\frac{6272}{6561}=1-\left(\frac{17}{81}\right)^2.
\]

It also evaluates the actual square-family Bellman recursion, using exact
factored fibers rather than assuming tensorization, for both consecutive and
alternating 14-coordinate orders. Each 256-bit Arb result agrees with twice
the corresponding base-order optimum; the residual enclosures have radius
below \(1.4\times10^{-70}\). The alternating check visits 2,985,831 states and
hulls 405,641 ambiguous clamps. These finite checks adversarially test the
load-bearing identity. The all-order and all-\(k\) conclusion follows from the
exact Bellman induction below, not from numerical extrapolation.

The factored representation used there still encodes the block structure whose
consequence is being tested. A separate audit therefore materializes the same
square as a plain list of 2,025 fourteen-bit masks and runs the generic
prefix-fiber evaluators on it, with no product-aware shortcut
([`audit_square_direct.py`](audit_square_direct.py)). For five declared global
orders — consecutive, block-swapped, reversed, alternating, and a seeded
interleaving — the directly computed fixed-order \(C_{+,\pi}\) and \(Q_\pi\)
equal the induced block sums to float64 rounding, worst gaps
\(8.9\times10^{-16}\) and \(3.6\times10^{-15}\), while the values themselves
take three distinct magnitudes, so the check tracks real order dependence. A
further resumable search
([`audit_tensorization.py`](audit_tensorization.py)) attacks additivity on 292
ordered products, 205 of them with unequal block dimensions and including dead
columns, duplicate columns, singletons, full power sets, and seeded random
pairs, over every global order of each product: no counterexample, worst gap
\(1.8\times10^{-15}\), and zero failures of the exact defect product.

That evidence was float64, so it could not separate exact additivity from
additivity up to machine epsilon. The exact rational evaluator closes that gap:
[`audit_tensorization_exact.py`](audit_tensorization_exact.py) repeats the
attack with certified rational enclosures on 35 ordered products — 20 with
unequal block dimensions, including singletons, full power sets, dead columns,
duplicate columns, deterministic coordinates, union-closed against non
union-closed factors, and seeded random pairs — over **every** global order of
every product, 15,010 orders in total. For each order it demands that the
enclosure of \(C_{+,\pi}(\mathcal H)\) intersect the enclosure of
\(C_{+,\pi_F}(\mathcal F)+C_{+,\pi_G}(\mathcal G)\); non-intersection would be
a counterexample rather than rounding. The result is zero intersection
failures, zero failures of the exact defect identity, worst \(Q\) discrepancy
\(1.4\times10^{-27}\) — the certificate width — and a worst fixed-order Bellman
discrepancy of **exactly zero**.

The computational lemma supplies a strictly negative base and direct product
controls. The passage to every \(k\) is the mathematical product argument
below.

### Independent clamp-free certificate

The strict negativity of the base has a second certificate that does not use
Arb, the original Bellman implementation, or any clamp classification. For a
fixed order, define \(W_\pi\) by the same terminal condition as
\(C_{+,\pi}\). At a state where all four children occur, replace the feasible
action interval by all of \([0,1]\). If the child values are \(W_{ac}\), this
relaxed step is

\[
W=-W_{10}r-W_{01}p+W_{00}+W_{11}(p+r)
  +\log_2\!\left(1+2^D\right),
\qquad
D=-W_{00}+W_{10}+W_{01}-W_{11}.
\tag{R}
\]

At a degenerate state the coupling action is forced, so the ordinary boundary
recurrence is used. Backward induction gives
\(C_{+,\pi}(\mathcal B)\le W_\pi(\mathcal B)\). Indeed, before enlarging the
action interval, every feasible transition coefficient is nonnegative and the
four coefficients sum to one. Replacing each true child by an upper child
therefore preserves the inequality. Maximizing the resulting concave function
over the larger interval \([0,1]\) can only increase it, and its maximum is
exactly (R).

The standalone standard-library checker
[`verify_gate_b_dyadic.py`](verify_gate_b_dyadic.py) evaluates \(Q\) and the
average of \(W_\pi\) with 160-bit outward-rounded dyadic intervals. Its
logarithm uses the rational atanh series with an exact tail bound; its
exponential uses range reduction, a rational Taylor tail bound, and repeated
squaring. It reconstructs the family from the seven cells and exactly verifies
the 240-element automorphism group and its 21 equal order orbits. It obtains

\[
\begin{aligned}
Q(\mathcal B)
&<5.447354977380321988160699884419,\\
\mathbb E_\pi W_\pi(\mathcal B)
&<5.933885055704082998831935023544,\\
A_+(\mathcal B)
&<-0.0271742911034863966862278304907
<-\frac1{40}.
\end{aligned}
\]

Thus this arithmetically independent, deliberately looser certificate already
proves the theorem: its powers have repair ratio

\[
\frac{-A_+(\mathcal F_k)}{\varepsilon_\vee(\mathcal F_k)}
>\frac{k/40}{1-(17/81)^k}>\frac{k}{40}\longrightarrow\infty.
\]

The stronger Arb bound \(-7/250\) is retained for the sharp finite-base record,
but the unboundedness conclusion no longer depends on its clamp decisions or
on the `python-flint` arithmetic implementation.

### Independent second base without an order quotient

The conclusion also follows from a different family, so it does not depend on
the seven-coordinate candidate being the right reconstruction of the search
output.  Partition six coordinates into two triples and let \(\mathcal D\) be
the union of the five weight cells

\[
(0,0),(0,1),(1,0),(1,2),(2,1).
\]

The third-party-free checker
[`verify_gate_b_n6_dyadic.py`](verify_gate_b_n6_dyadic.py) reconstructs all 25
rows from this definition.  It verifies exactly that every coordinate count is
10, the total incidence is 60, the exact Reimer threshold is 59, the columns
are active and distinct, and 444 of the \(25^2\) ordered joins are missing.
Thus

\[
\varepsilon_\vee(\mathcal D)=\frac{444}{625},
\qquad
1-\varepsilon_\vee(\mathcal D)=\frac{181}{625}.
\]

Unlike the seven-coordinate certificate, this computation evaluates all
\(6!=720\) coordinate orders and uses no automorphism quotient.  Applying the
same clamp-free relaxation as above with 160-bit outward-rounded dyadic
intervals gives

\[
\begin{aligned}
Q(\mathcal D)
&<4.616383859857110260038923460896,\\
\mathbb E_\pi W_\pi(\mathcal D)
&<5.034353951629763221689815640761,\\
A_+(\mathcal D)
&<-0.0125897106568747589615082453240
<-\frac1{80}.
\end{aligned}
\]

This route shares the already-audited dyadic transcendental primitive with the
preceding certificate, but it shares neither the finite family nor an order
symmetry reduction.  A second checker,
[`verify_gate_b_n6_arb.py`](verify_gate_b_n6_arb.py), passes the same family
through the generic 256-bit Arb Bellman evaluator and its exact 72-element
automorphism group.  It independently certifies the stronger bound
\[
A_+(\mathcal D)<-\frac{17}{1250}.
\]
Thus the second construction has both an all-order clamp-free certificate and
an arithmetically disjoint clamp-classified certificate.

For \(\mathcal G_k=\mathcal D^{\boxtimes k}\), the product lemmas below give

\[
A_+(\mathcal G_k)=kA_+(\mathcal D),\qquad
\varepsilon_\vee(\mathcal G_k)
=1-\left(\frac{181}{625}\right)^k.
\]

Every coordinate occurs in \(10\cdot25^{k-1}=(2/5)25^k\) rows.  The average row
size is \(12k/5\), while \(25^5<2^{24}\) implies
\(\frac12\log_2(25^k)<12k/5\), so every power is cap/Reimer and normalized.
Consequently this second construction alone proves

\[
\frac{-A_+(\mathcal G_k)}{\varepsilon_\vee(\mathcal G_k)}
>
\frac{k/80}{1-(181/625)^k}
>\frac{k}{80}\longrightarrow+\infty.
\]

### Exact rational certificate with no relaxation and no interval library

Both preceding routes bound \(A_+\) from above, but each asks for extra trust:
the Arb route uses a third-party arithmetic library and a three-way clamp
classification, and the dyadic route enlarges every feasible action interval to
\([0,1]\), so it certifies a relaxation rather than \(C_+\) itself.

[`verify_gate_b_rational.py`](verify_gate_b_rational.py) removes both
concessions. It evaluates the recurrence of
[DEFINITIONS.md](DEFINITIONS.md) over the **exact** feasible interval
\([s^*(p,r),U(p,r)]\) in ordinary rational arithmetic, and it produces a
*two-sided* enclosure. It needs exactly three certified scalar primitives.

**Primitive 1 (binary logarithm).** For \(y>0\) and any bit \(b\in\{0,1\}\),

\[
\log_2y=\frac{b+\log_2\!\left(y^2/2^{b}\right)}{2}.
\]

Iterating \(N\) times and using monotonicity of \(\log_2\), replacing each
intermediate state by a rational upper (lower) bound keeps the accumulated
value an upper (lower) bound. Choosing \(b=1\) exactly when the rounded square
is at least two keeps every state in \((0,2]\) upward and in \([1,2]\)
downward, so the discarded tail \(2^{-N}\log_2y_N\) lies in
\([0,2^{-N}]\). Adding \(2^{-N}\) upward and nothing downward therefore gives
certified rational bounds. All arithmetic is integer.

**Primitive 2 (binary exponential).** With \(r_0=2\) and
\(r_{j+1}=\lceil\sqrt{r_j}\rceil\) rounded up to a fixed dyadic denominator,
induction gives \(r_j\ge2^{2^{-j}}\). Rounding a rational exponent up to \(N\)
binary places and multiplying the \(r_j\) over its set bits therefore
upper-bounds \(2^x\). Only `math.isqrt` is used.

**Primitive 3 (entropy and its constrained maximum).** Since both weights in
\(h(s)=-s\log_2s-(1-s)\log_2(1-s)\) are negative, a lower bound on each
logarithm gives an upper bound on \(h\), and conversely. For the Bellman step,
\(h(s)+cs\) is concave with derivative \(\log_2((1-s)/s)+c\), so on
\([\ell,u]\subseteq[0,1]\):

* if \(c\le\log_2(\ell/(1-\ell))\) the maximum is exactly at \(\ell\);
* if \(c\ge\log_2(u/(1-u))\) it is exactly at \(u\);
* otherwise it is at most the unconstrained Fenchel value
  \(\max_{s\in[0,1]}[h(s)+cs]=\log_2(1+2^{c})\).

Each branch is a valid upper bound, and the endpoint branches are exact. The
two sign tests are decided with certified bounds of the same logarithm, and the
third branch is always available, so the case split never needs an exact
transcendental comparison.

**Monotone induction.** On the feasible interval every transition probability
\(P_{ac}(s)\) is nonnegative and \(\sum_{ac}P_{ac}(s)=1\), so a Bellman step is
monotone in its children. Replacing children by upper bounds, and then the
affine data \((\text{const},c)\) by upper bounds — legitimate because
\(s\ge0\) — keeps the result an upper bound of the true optimum. For the lower
bound the same recursion evaluates one explicitly *feasible* rational action
near \(1/(1+2^{-c})\) against pessimistic children; that is the value of an
admissible one-sided coupling, hence at most \(C_{+,\pi}\).

**Result.** Evaluating all \(6!=720\) orders of \(\mathcal D\) and all
\(7!=5040\) orders of \(\mathcal B\) — no automorphism quotient in either case —
gives the certified rational enclosures

\[
\begin{aligned}
A_+(\mathcal D)&\in
[-0.0136721077321777735618599156610,\,
  -0.0136721077321777735618599129978],\\
A_+(\mathcal B)&\in
[-0.0286491867944683178160269819327,\,
  -0.0286491867944683178160269764548].
\end{aligned}
\]

Both enclosures have width below \(6\times10^{-27}\) and both lie strictly
below \(-17/1250\) and \(-7/250\) respectively, so this single artifact
re-derives the sharp form of *both* base lemmas. The Arb values quoted earlier
lie strictly inside these enclosures, which also confirms that retaining the
unconstrained maximum at the 982 and 305 clamp-ambiguous states cost less than
\(10^{-26}\).

The upper endpoints are what the theorem consumes. The lower endpoints are not
needed for unboundedness; they matter because a sign error, an inverted clamp
comparison, or a mis-stated feasible interval would break the sandwich instead
of silently shifting one bound.


## Product lemmas

Let \(\mathcal F\subseteq2^{[n]}\) and
\(\mathcal G\subseteq2^{[d]}\), with uniform laws, and put
\(\mathcal H=\mathcal F\boxtimes\mathcal G\).

### Lemma 1: the closure defect multiplies through success

A joined product row belongs to \(\mathcal H\) if and only if its
\(\mathcal F\)-block join belongs to \(\mathcal F\) and its
\(\mathcal G\)-block join belongs to \(\mathcal G\). The two block events are
independent under two independent uniform product rows. Hence

\[
1-\varepsilon_\vee(\mathcal H)
=(1-\varepsilon_\vee(\mathcal F))
 (1-\varepsilon_\vee(\mathcal G)).
\tag{1}
\]

Iteration gives the displayed formula for \(\mathcal F_k\).

### Lemma 2: the Shapley iid cost is additive

Fix a global coordinate order \(\pi\) and let \(\pi_F,\pi_G\) be its induced
relative orders on the two blocks. For a coordinate in the \(\mathcal F\)
block, conditioning on earlier \(\mathcal G\) coordinates does not change the
conditional-one probability in \(\mathcal F\), because a uniform product row
has independent blocks. Its local iid OR-entropy is therefore exactly the
corresponding local term for \(\mathcal F\) under \(\pi_F\). The same argument
holds in the other block, so

\[
Q_\pi(\mathcal H)=Q_{\pi_F}(\mathcal F)+Q_{\pi_G}(\mathcal G).
\tag{2}
\]

Under a uniform global order, each induced relative order is uniform. Averaging
(2) gives

\[
Q(\mathcal H)=Q(\mathcal F)+Q(\mathcal G).
\tag{3}
\]

### Lemma 3: the one-sided Bellman cost is additive

Again fix \(\pi\). At any pair of product-row prefixes, split the state into
its \(\mathcal F\) and \(\mathcal G\) components. We prove by backward
induction on the remaining global suffix that its Bellman value is

\[
V_H=V_F+V_G,
\tag{4}
\]

where each factor value uses the induced remaining suffix and its own two
prefixes.

The identity is zero equals zero at the terminal state. Suppose the next
coordinate lies in the \(\mathcal F\) block. Its conditional marginals
\(p,r\), feasible interval \([s^*(p,r),U(p,r)]\), and four transition
probabilities depend only on the two \(\mathcal F\) prefixes. By the induction
hypothesis, every child continuation value is its \(\mathcal F\)-child value
plus the same unchanged number \(V_G\). Since the four transition
probabilities sum to one, \(V_G\) leaves the maximization as an additive
constant. The remaining maximization is exactly the \(\mathcal F\) Bellman
step. This proves (4); the \(\mathcal G\) case is identical.

Explicitly, for every feasible action \(s\),

\[
\sum_{a,c}P_{ac}(s)
=(1-s)+(s-r)+(s-p)+(p+r-s)=1,
\]

and therefore

\[
\max_s\!\left[h(s)+\sum_{a,c}P_{ac}(s)
  (V_{F,ac}+V_G)\right]
=V_G+\max_s\!\left[h(s)+\sum_{a,c}P_{ac}(s)V_{F,ac}\right].
\]

The common shift also cancels from the optimizer's slope:

\[
-(V_{F,00}+V_G)+(V_{F,10}+V_G)+(V_{F,01}+V_G)
 -(V_{F,11}+V_G)
=-V_{F,00}+V_{F,10}+V_{F,01}-V_{F,11}.
\]

In particular, allowing a policy at one block to observe prefixes from the
other block cannot increase the value: the other-block term cancels from the
Bellman slope and action choice. At the root,

\[
C_{+,\pi}(\mathcal H)
=C_{+,\pi_F}(\mathcal F)+C_{+,\pi_G}(\mathcal G).
\tag{5}
\]

Averaging over global orders gives

\[
C_+(\mathcal H)=C_+(\mathcal F)+C_+(\mathcal G).
\tag{6}
\]

### Corollary: \(A_+\) is additive

Since \(|\mathcal H|=|\mathcal F||\mathcal G|\), equations (3), (6), and
\(\log_2|\mathcal H|=\log_2|\mathcal F|+\log_2|\mathcal G|\) imply

\[
A_+(\mathcal F\boxtimes\mathcal G)
=A_+(\mathcal F)+A_+(\mathcal G).
\tag{7}
\]

This identity uses the repository's fixed \(\alpha\) without changing its
value.

## Admissibility of every power

For \(\mathcal F_k=\mathcal B^{\boxtimes k}\):

- \(|\mathcal F_k|=45^k\) and the dimension is \(7k\).
- A coordinate belongs to
  \(18\cdot45^{k-1}=(2/5)45^k\) rows, exactly meeting the cap.
- The total incidence is \(126k45^{k-1}\), so the average row size is
  \(14k/5\). The integer inequality
  \(45^5=184{,}528{,}125<268{,}435{,}456=2^{28}\) gives
  \(\log_2 45<28/5\), and therefore

  \[
  \frac{14k}{5}>\frac12\log_2(45^k).
  \]

  Reimer holds strictly.
- Activity and separation within each block are inherited from \(\mathcal B\).
  Columns in different blocks are distinct because every base column takes
  both values and product blocks vary independently. Hence every power is
  normalized.

Powers beyond the square were previously only integer arithmetic derived from
these formulas. [`audit_power_admissibility.py`](audit_power_admissibility.py)
now instantiates \(\mathcal F_k\) for \(k\le3\) and recomputes every claim from
the constructed rows, counting missing ordered joins by brute force. At
\(k=3\): 91,125 rows, all 21 coordinate counts equal to the cap bound 36,450,
incidence 765,450 against Reimer threshold 750,668, active and separating, and
8,227,000,000 of 8,303,765,625 ordered joins missing, i.e. exactly
\(526528/531441=1-(17/81)^3\). That count uses no product shortcut, so it is an
independent confirmation of Lemma 1 one power beyond the square.

Combining admissibility, (1), (7), and the certified base lemma completes the
proof.
Every constructed defect lies strictly between zero and one, so the proof
never divides by zero. The result is an all-\(k\) consequence of the product
lemmas; no finite-\(n\) optimality or extrapolation is used.

The only numerical input is the sign of \(A_+(\mathcal B)\), and it is now
established three times by arithmetically disjoint routes: the 256-bit Arb
certificate with its clamp classification, the standard-library dyadic
clamp-free relaxation, and the exact rational two-sided enclosure over the true
feasible action set. The algebra those routes share — unit transition mass,
shift invariance of the Bellman slope, and the closed form
\(\max_s[sD+h(s)]=\log_2(1+2^D)\) with strictly negative second derivative —
is machine-checked symbolically in `test_gate_b.py`.

All three certificates are ordinary verified computation, not proof-assistant
artifacts. The product lemmas are human-audited mathematics supported by the
falsification searches above.

## Corollary: the divergence is exactly linear in the dimension

Unboundedness alone does not say how fast the obstruction grows. Define

\[
c_{\rm cl}^\star(n)=
\sup\Big\{\tfrac{-A_+(\mathcal F)}{\varepsilon_\vee(\mathcal F)}:
\mathcal F\subseteq2^{[n']},\ n'\le n,\ \mathcal F\text{ cap/Reimer},\
A_+(\mathcal F)<0\Big\}.
\]

**Lower bound.** For a certified base of dimension \(d\) with
\(A_+\le-\delta<0\), the powers \(\mathcal F_k\) for \(k=\lfloor n/d\rfloor\)
live in \(2^{[dk]}\subseteq2^{[n]}\) and give

\[
c_{\rm cl}^\star(n)\ \ge\ \frac{\delta\lfloor n/d\rfloor}
{1-(1-\varepsilon_\vee)^{\lfloor n/d\rfloor}}
\ >\ \delta\left\lfloor\frac nd\right\rfloor .
\]

With \(\mathcal B\) this is \(c_{\rm cl}^\star(n)>\frac7{250}\lfloor n/7\rfloor\),
i.e. slope \(1/250\) per coordinate; with \(\mathcal D\) it is
\(\frac{17}{1250}\lfloor n/6\rfloor\), slope \(17/7500\). Both are exact
rational consequences of the enclosures above.

**Matching upper bound under a defect floor.** Every local term of \(Q\) is a
binary entropy, so \(Q\ge0\); \(C_+\) is a maximum over policies of a sum of
binary entropies, so \(C_+\ge0\). With \(0<\alpha<1\) this gives

\[
A_+(\mathcal F)=(1-\alpha)Q+\alpha C_+-\log_2m\ \ge\ -\log_2m\ \ge\ -n ,
\]

because \(m=|\mathcal F|\le2^{n}\). Hence for any admissible family on at most
\(n\) coordinates whose defect satisfies \(\varepsilon_\vee\ge\varepsilon_0>0\),

\[
\frac{-A_+(\mathcal F)}{\varepsilon_\vee(\mathcal F)}
\ \le\ \frac{\log_2m}{\varepsilon_0}\ \le\ \frac n{\varepsilon_0}.
\]

So on the subclass with defect bounded away from zero the Gate B ratio grows
**exactly of order \(n\)**: the two certified families realize slopes
\(1/250\) and \(17/7500\), and no family can exceed \(1/\varepsilon_0\) per
coordinate. In particular the divergence proved above is not an artifact of a
wild sequence; it is the true growth rate up to a constant factor, and the only
way to escape it is to force \(\varepsilon_\vee\to0\), which is precisely the
separate local-stability question left open below.

A remark on the unrestricted supremum at fixed \(n\): dropping the defect floor
leaves only \(\varepsilon_\vee\ge1/m^2\) for a non-union-closed family, so the
trivial bound degrades to \(n\,m^2\le n4^{n}\). The theorem does not claim that
this is attained; the linear statement above is the sharp one available with
these definitions.

## Scope and stronger interpretation

This theorem resolves the repository's displayed universal supremum. The
constructed defects tend to one, not zero. Therefore it does **not** resolve a
different local-stability question restricted to sequences with
\(\varepsilon_\vee\to0\). The prior plan's phrase "near-UC adversarial search"
was a method proposal; no such restriction appeared in the Gate B definition.
The distinction is recorded explicitly rather than silently changing Gate B.

The result also does not produce a union-closed counterexample or prove
Frankl's conjecture. It proves that every fixed scalar multiple of the ordered
pair closure defect is structurally incapable of repairing this particular
\(A_+\) relaxation over all cap/Reimer families.
