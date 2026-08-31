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
leaves only Lemma 4 below, \(\varepsilon_\vee\ge2/m^2\), for a
non-union-closed family, so the trivial bound degrades to
\(nm^2/2\le n4^{n}/2\). The theorem does not claim that this is attained; the
linear statement above is the sharp one available with these definitions.

## The local regime

The corollary just proved isolates one escape and only one: forcing
\(\varepsilon_\vee\to0\). This section makes that regime precise, proves what is
elementary about it, and states exactly where it is blocked. For
\(0<\varepsilon\le1\) put

\[
c_{\rm loc}(\varepsilon)=
\sup\Big\{\tfrac{-A_+(\mathcal F)}{\varepsilon_\vee(\mathcal F)}:
\mathcal F\text{ cap/Reimer},\ A_+(\mathcal F)<0,\
0<\varepsilon_\vee(\mathcal F)\le\varepsilon\Big\},
\qquad
c_{\rm loc}=\lim_{\varepsilon\to0^+}c_{\rm loc}(\varepsilon),
\]

with \(\sup\emptyset=0\). The limit exists in \([0,+\infty]\) because
\(c_{\rm loc}\) is nondecreasing in \(\varepsilon\). Gate B is the statement
\(c_{\rm loc}(1)=+\infty\); the open question is the value of \(c_{\rm loc}\).

### Lemma 4: failures come in pairs

If \(\mathcal F\) is not union-closed then
\(\varepsilon_\vee(\mathcal F)\ge2/m^2\).

*Proof.* Let \(S=\{(X,Y)\in\mathcal F^2:X\cup Y\notin\mathcal F\}\), which is
nonempty by assumption. \(S\) is invariant under exchanging the two
coordinates, and it contains no diagonal element because \(X\cup X=X\in\mathcal
F\). Hence \(S\) is a disjoint union of two-element orbits, so \(|S|\ge2\) is
even. \(\square\)

### Lemma 5: products never dilute the defect

For \(\mathcal H=\mathcal F\boxtimes\mathcal G\),

\[
1-\varepsilon_\vee(\mathcal H)
=\big(1-\varepsilon_\vee(\mathcal F)\big)\big(1-\varepsilon_\vee(\mathcal G)\big),
\qquad\text{hence}\qquad
\varepsilon_\vee(\mathcal H)\ \ge\
\max\{\varepsilon_\vee(\mathcal F),\varepsilon_\vee(\mathcal G)\},
\]

with equality in the second display if and only if the other factor is
union-closed.

*Proof.* The identity is Lemma 1. Monotonicity follows since
\(1-\varepsilon_\vee(\mathcal H)=(1-\varepsilon_\vee(\mathcal F))(1-\varepsilon_
\vee(\mathcal G))\le1-\varepsilon_\vee(\mathcal F)\), and equality forces
\(\varepsilon_\vee(\mathcal G)=0\). \(\square\)

**Consequence for the proof above.** Every family in the Gate B construction
satisfies \(\varepsilon_\vee(\mathcal B^{\boxtimes k})\ge64/81\) and
\(\varepsilon_\vee(\mathcal D^{\boxtimes k})\ge444/625\). So the Cartesian-power
mechanism cannot exhibit even one family with
\(\varepsilon_\vee<444/625=0.7104\), and therefore contributes nothing to
\(c_{\rm loc}\). The local question is untouched by the theorem, not merely
unaddressed by its statement.

### Lemma 6: the local regime forces large families

For a non-union-closed admissible \(\mathcal F\),

\[
\frac{-A_+(\mathcal F)}{\varepsilon_\vee(\mathcal F)}
\ \le\ \frac{m^2\log_2m}{2},
\]

and \(\varepsilon_\vee(\mathcal F)\le\varepsilon\) forces
\(m\ge\sqrt{2/\varepsilon}\).

*Proof.* Combine \(-A_+\le\log_2m\), proved in the corollary above from
\(Q,C_+\ge0\), with Lemma 4. The second claim is Lemma 4 rearranged.
\(\square\)

So \(c_{\rm loc}(\varepsilon)\) is finite on every subclass of bounded size, and
any sequence witnessing \(c_{\rm loc}=+\infty\) must have \(m\to\infty\), hence
\(n\to\infty\) as well since \(m\le2^n\).

### Proposition 7: the zero-defect boundary is Frankl's conjecture at 2/5

An admissible family with \(\varepsilon_\vee=0\) is union-closed and has every
degree at most \(\lfloor2m/5\rfloor<m/2\), i.e. it is a counterexample to
Frankl's union-closed sets conjecture. Two consequences.

*(i)* **[REPORTED]** Frankl's conjecture is verified for \(n\le12\)
(Vučković and Živković, *The 12-Element Case of Frankl's Conjecture*, IPSI BgD
Transactions on Internet Research 13(1), 65-71, 2017), which together with
Faro's bound and the Roberts-Simpson estimate \(|\mathcal F|\ge4q-1\) for a
minimal counterexample on \(q\) elements gives the conjecture for \(m\le50\)
(as stated in arXiv:1711.04276). Neither range is re-verified here.
Hence no admissible family on at most twelve coordinates has zero defect: every
search reported below runs at \(n\le8\), strictly inside the verified range, so
its families all have \(\varepsilon_\vee\ge2/m^2>0\) by Lemma 4 and the
extended-value convention never activates.

*(ii)* Suppose one could prove a universal positive defect floor,

\[
\varepsilon_\vee(\mathcal F)\ \ge\ \varepsilon_*>0
\qquad\text{for every cap-}2/5\text{ family }\mathcal F .
\]

Then no union-closed family satisfies the cap, i.e. every union-closed family
has an element in more than \(2m/5\) of its sets. The best published constant
for that statement is \(\psi=(3-\sqrt5)/2=0.381966\ldots\), with small explicit
improvements above it; this repository's own Campaign I certified
\(\psi+10^{-4}\). A positive floor at cap \(2/5\) would therefore improve the
union-closed frontier by about \(0.018\), which is far beyond current technique.

**This is the precise blocker.** Emptying the local regime by a theorem is at
least as hard as a major advance on Frankl's conjecture. Populating it instead
requires exhibiting an admissible family that is simultaneously nearly union
closed and has \(A_+<0\); by Lemma 5 no product of the certified bases is such a
family, and by Lemma 6 any such family must be large. Both directions are
recorded as open, with the searched evidence in
[EXPERIMENTS.md](EXPERIMENTS.md).

### Proposition 8: the negativity region reaches defect 1336/2025

Let \(\mathcal L\) be the 45-row family on seven coordinates listed in
[`certificates/gate_b_lowdefect_rational_v1.json`](certificates/gate_b_lowdefect_rational_v1.json).
Then \(\mathcal L\) is normalized and cap/Reimer-admissible, every coordinate
degree equals the cap \(18\), the incidence is \(126\ge124\), and

\[
\varepsilon_\vee(\mathcal L)=\frac{1336}{2025},
\qquad
A_+(\mathcal L)\in
[-0.00085518533372314171100789579,\,-0.00085518533372314171100789029],
\]

so \(A_+(\mathcal L)<-1/1200<0\). Consequently

\[
e^\star:=\min\{\varepsilon_\vee(\mathcal F):
\mathcal F\text{ admissible},\ A_+(\mathcal F)<0\}
\ \le\ \frac{1336}{2025}=0.6597\ldots
\]

*Proof.* The combinatorial facts are exact integer computations on the row
list, replayed by the checker before it evaluates anything. The enclosure is
produced by the exact rational evaluator of the previous section, applied to
this family over **all** \(5040\) coordinate orders; its automorphism group is
trivial, so the run records 5040 distinct enclosures and no symmetry quotient is
used or available. The upper endpoint is a rational number smaller than
\(-1/1200\). \(\square\)

By Lemmas 2, 3 and 5 the powers \(\mathcal L^{\boxtimes k}\) are admissible with
\(A_+=kA_+(\mathcal L)\) and \(\varepsilon_\vee=1-(689/2025)^k\), so they give a
third infinite construction with ratio exceeding \(k/1200\). Its role is not
sharpness -- both published bases are far more negative -- but *defect*: it
certifies that negativity survives down to closure defect \(0.6598\), against
\(64/81=0.7901\) for \(\mathcal B\) and \(444/625=0.7104\) for \(\mathcal D\).

Together with the searched evidence in [EXPERIMENTS.md](EXPERIMENTS.md), which
found no negative family below defect \(0.64\) at \(m=45\) and a least
\(A_+\) rising monotonically to \(+0.0795\) as the defect cap falls to
\(0.35\), the current picture of the local regime is an interval of large defect
on which \(A_+<0\), with \(-A_+\to0\) at its lower endpoint. That is the
behaviour of a *bounded* local ratio. The certified frozen target gives
\(-A_+(\mathcal L)/\varepsilon_\vee(\mathcal L)>27/21376=0.001263\ldots\) as a
lower bound, and the enclosure endpoints pin the value itself to
\(0.0012962\ldots\), twenty-eight times below the \(0.0362591\ldots\) of
\(\mathcal B\).
Nothing here decides \(c_{\rm loc}\), and by Proposition 7 deciding it downward
is at least as hard as a major advance on Frankl's conjecture.

## The size ceiling, the defect floor, and the eighth coordinate

### Lemma 9: the admissible size ceiling

Every admissible \(\mathcal F\subseteq2^{[n]}\) satisfies

\[
\log_2 m\ \le\ \frac{4n}{5},
\qquad\text{equivalently the integer test}\qquad m^5\le 2^{4n}.
\]

*Proof.* Incidence counted by coordinates is
\(\sum_{A\in\mathcal F}|A|=\sum_{i\in[n]}\deg_i\), and the cap bounds every
degree by \(\lfloor 2m/5\rfloor\le 2m/5\), so the incidence is at most
\(2nm/5\). Reimer admissibility demands incidence at least \(m\log_2m/2\).
Hence \(m\log_2m/2\le 2nm/5\); divide by \(m>0\). \(\square\)

The bound is close to sharp: the largest admissible sizes are \(45\) at
\(n=7\) (\(\log_2 45=5.4919\) against \(5.6\)) and \(445\) at \(n=11\)
(\(8.7977\) against \(8.8\)). The repository already carried this inequality as
the per-base `reimer_witness` string; Lemma 9 is the statement that it holds for
*every* admissible family and depends only on \(n\).

### Corollary 10: the numerator grows exactly linearly

\(Q\ge0\) and \(C_+\ge0\) give \(-A_+=\log_2m-(1-\alpha)Q-\alpha C_+\le\log_2m\),
so by Lemma 9 every admissible family on \(n\) coordinates has

\[
-A_+(\mathcal F)\ \le\ \frac{4n}{5}.
\]

The certified powers give \(-A_+(\mathcal B^{\boxtimes k})=k\,(-A_+(\mathcal B))\)
on \(7k\) coordinates. Writing \(\Lambda(n)\) for the supremum of \(-A_+\) over
admissible families on at most \(n\) coordinates,

\[
\frac{7}{250}\left\lfloor\frac n7\right\rfloor\ \le\ \Lambda(n)\ \le\ \frac{4n}{5},
\]

so \(\Lambda(n)=\Theta(n)\), bounded on both sides. Only the denominator of the
repair ratio is still open.

### Lemma 11: the union-growth defect floor

Let \(M=\max_{A\in\mathcal F}|A|<n\) and \(\bar s=\frac1m\sum_{A}|A|\). Then

\[
\varepsilon_\vee(\mathcal F)\ \ge\ \frac{\frac85\bar s-M}{\,n-M\,}.
\]

*Proof.* With \(p_i=\deg_i/m\) the cap gives \(p_i\le2/5\), so
\[
\mathbb E|X\vee Y|=\sum_i\bigl(1-(1-p_i)^2\bigr)=\sum_i p_i(2-p_i)
\ \ge\ \sum_i p_i\left(2-\tfrac25\right)=\tfrac85\bar s .
\]
A union that lands in \(\mathcal F\) has size at most \(M\), and any union has
size at most \(n\), so \(\mathbb E|X\vee Y|\le(1-\varepsilon_\vee)M+\varepsilon_\vee n\).
Combine and solve for \(\varepsilon_\vee\). \(\square\)

Reimer gives \(\bar s\ge\frac{\log_2m}{2}\) and hence \(\frac85\bar s\ge\frac45\log_2m\).

### Corollary 12: without a dominant set the ratio is O(n)

If \(M\le\theta\cdot\frac85\bar s\) for some \(\theta<1\), then by Corollary 10
and Lemma 11

\[
\frac{-A_+}{\varepsilon_\vee}\ \le\ \frac{(4n/5)\,(n-M)}{\frac85\bar s-M}
\ \le\ \frac{5n}{4(1-\theta)} .
\]

So the repair ratio is \(O(n)\) on every family without a dominant set, matching
the certified \(\Theta(n)\) lower bound. **Any sequence whose ratio diverges
faster than linearly, and in particular any approach to the local regime, must
contain sets of size at least \(\frac85\bar s\ge\frac45\log_2m\).** The
\(n=6\) base has \(M=3\) below its threshold \(96/25=3.84\), and its floor
\(7/25\) is active; both \(n=7\) bases have \(M\in\{6,7\}\) above their threshold
\(112/25=4.48\), and every \(n=8\) witness contains \([8]\) itself. The
dominant set is not incidental, it is what the low-defect families are made of.

### Lemma 13: the downset capacity floor

With \(N(A)=|\{B\in\mathcal F:B\subseteq A\}|\),

\[
\varepsilon_\vee(\mathcal F)\ \ge\ 1-\frac1{m^2}\sum_{A\in\mathcal F}N(A)^2 .
\]

*Proof.* \(X\vee Y=A\) forces \(X,Y\in\mathcal F\) with \(X,Y\subseteq A\), so at
most \(N(A)^2\) ordered pairs have union \(A\); sum over \(A\in\mathcal F\) and
compare with \((1-\varepsilon_\vee)m^2\). \(\square\)

This one is attained: the exhaustive audit over all 366 admissible families with
\(n\le4\) reports tightest slack exactly \(0\).

### Proposition 14: the certified bases admit no addition at all

At \(n=6\) and \(n=7\) the largest admissible sizes are \(25\) and \(45\), and
\(\mathcal D\), \(\mathcal B\), \(\mathcal L\) all attain them with *every*
coordinate degree exactly at the cap. Consequently, for every
\(C\subseteq[n]\) not already present, \(\mathcal F\cup\{C\}\) is inadmissible:
a nonempty \(C\) pushes some degree past \(\lfloor2(m+1)/5\rfloor\), which does
not grow at \(m=25\) or \(m=45\), and \(C=\emptyset\) leaves the incidence fixed
while \(R_{m+1}\) rises above it.

*Verified exhaustively* over all \(2^n\) candidate sets for all three bases by
[`bound_local_regime.py`](bound_local_regime.py), which reports zero admissible
additions. So the "add the missing unions" route to lower defect -- the one
amplification mechanism that lowers \(\varepsilon_\vee\) and raises \(\log_2 m\)
simultaneously -- is not merely unpromising at the certified bases, it is empty.
Lowering the defect requires a larger ground set, which is exactly where the
next section goes.

### Proposition 15: the local frontier at eight coordinates

The complete \(S_2\times S_6\) block-symmetric class at \(n=8\) contains 2,272
canonical admissible families, of which 21 have negative float64 objective and
**20 are certified negative in exact rational arithmetic**. Among them:

* a separating 75-row family with \(\varepsilon_\vee=1144/1875=0.610133\ldots\)
  and \(A_+\le-0.002032376\), improving the certified local frontier from
  \(1336/2025=0.659753\ldots\);
* a separating 75-row family with \(\varepsilon_\vee=468/625\) and
  \(A_+\le-0.027744001\), whose repair ratio \(0.0370512\ldots\) is the first in
  this repository to exceed the published base's \(0.0362591\ldots\);
* two non-separating 70-row families reaching \(\varepsilon_\vee=139/245=0.567347\ldots\)
  and ratio \(0.0580391\ldots\). `DEFINITIONS.md` records normalization as a
  search and reporting condition rather than a condition in the displayed
  supremum, so these are admissible; the flag is carried with every claim.

Size 70 and 75 are **infeasible at \(n=7\)**: \(7\lfloor2\cdot70/5\rfloor=196\)
is below \(R_{70}=215\). The eighth coordinate supplies the missing incidence,
and in the two 70-row families it is an exact duplicate of the first. Coordinate
cloning, which cannot change the defect or the size, is therefore what makes
these sizes admissible at all: it is the amplification mechanism that Lemma 9
predicts, since raising \(n\) is the only way to raise the size ceiling.

**The census objective is a screen, not an estimate.** The 80-row family at
\(\varepsilon_\vee=2553/3200\) reads \(-0.00026687\) in the census and
\(+0.000454098\) exactly. This is *not* float64 imprecision: on the same order
set the float evaluator agrees with the exact one to \(10^{-14}\), and the exact
orbit route agrees with all \(8!\) orders exactly. The census averages
\(C_{+,\pi}\) over one order per \(S_k\times S_{8-k}\) pattern, which represents
the objective only when the family's automorphism group *is* that block group.
For this class it is not: \(|\mathrm{Aut}|=1440\) matches \(|S_2\times S_6|\),
but only 120 of those elements preserve the block partition, so the 28
block-pattern orders cover just 12 of the 28 true orbits. Any subset average
lies within \(\alpha(\max_\pi C_{+,\pi}-\min_\pi C_{+,\pi})\) of the truth, at
most \(5.2\times10^{-3}\) over the registered bases, so candidates are re-ranked
exactly with a threshold of \(+0.01\); every number above comes from the exact
route.

### Proposition 16: no lemma above uses separation, and the growth constant improves

**No step of the product argument uses separation or activity.** Lemma 1 uses
only that the two block join events are independent under two independent
uniform product rows. Lemmas 2 and 3 use only that a uniform product row has
independent blocks, that the conditional marginals, feasible interval and four
transition probabilities at a coordinate depend on that block's prefixes alone,
and that the four transition probabilities sum to one. The Corollary uses only
\(|\mathcal H|=|\mathcal F||\mathcal G|\). In the admissibility section,
separation appears exactly once, in the last bullet, where it is *concluded* for
the powers of a normalized base; it is never a hypothesis.

Moreover admissibility of the powers needs nothing beyond the base's own
conditions. For a base with every degree equal to \(c\), incidence \(I\) and
size \(m\):

* *Cap.* The power has size \(m^k\) and coordinate degree \(cm^{k-1}\). Since
  \(c\le\lfloor2m/5\rfloor\le2m/5\), the degree is an integer at most
  \(2m^k/5\), hence at most \(\lfloor2m^k/5\rfloor\).
* *Reimer.* The power's incidence is \(kIm^{k-1}\), an integer, and the
  requirement is \(kIm^{k-1}\ge\lceil km^k\log_2m/2\rceil\). An integer
  dominating a real number dominates its ceiling, so it suffices that
  \(kIm^{k-1}\ge km^k\log_2m/2\), which cancels to \(2I\ge m\log_2m\), i.e. to
  the integer inequality \(m^m\le2^{2I}\) --- and that is precisely the base
  condition \(I\ge R_m\), because \(R_m\) is the least \(r\) with
  \(2^{2r}\ge m^m\). **One base-level Reimer check certifies every power.**

Consequently the two cloned-coordinate witnesses of Proposition 15 carry valid
infinite constructions. Certified in
[`certificates/gate_b_n8_clone_rational_v1.json`](certificates/gate_b_n8_clone_rational_v1.json),
with the powers audited by direct instantiation in
[`audit_clone_power.py`](audit_clone_power.py) (verdict
`EVERY_POWER_ADMISSIBLE`: degrees equal the cap exactly at every power, Reimer
strict, and at \(k=2\) the 4,900-row product recomputed from its rows has
defect \(210171/240100=1-(173/490)^2\) exactly):

\[
\mathcal C=\text{the 70-row }n=8\text{ witness},\qquad
\varepsilon_\vee(\mathcal C)=\frac{317}{490},\qquad
A_+(\mathcal C)\le-\frac{3}{80}.
\]

Its powers give \(c_{\rm cl}^\star(n)\ge\frac3{80}\lfloor n/8\rfloor\), i.e. an
asymptotic slope \(3/640=0.0046875\), against \(1/250=0.004\) for the published
\(n=7\) base --- an improvement by a factor \(75/64=1.1719\ldots\). Combined
with Corollary 10 the numerator is now pinned between

\[
\frac3{640}\,n-\frac3{80}\ \le\ \Lambda(n)\ \le\ \frac45\,n .
\]

The mechanism deserves a remark, because it is not the mechanism one would
guess. Cloning a coordinate changes neither the size, nor the closure defect,
nor the join structure; it only adds incidence. What it buys is *feasibility*:
size 70 is impossible at seven coordinates, where \(7\lfloor2\cdot70/5\rfloor=196\)
falls short of \(R_{70}=215\). So the eighth coordinate is spent entirely on
clearing Reimer, and the reward is a base whose \(-A_+\) per coordinate exceeds
anything available at seven. The cost is separation, which the definitions do
not charge for.

### Proposition 17: block-symmetric coverage at eight coordinates is complete

\(S_k\times S_{8-k}\) pairs \(k\) with \(8-k\), so \(k\in\{1,2,3,4\}\) exhausts
every block-symmetric class at \(n=8\). All four are enumerated and every
candidate whose screen value is below \(+0.01\) is exactly re-ranked:

| class | cells | raw masks | canonical families | screen-negative | certified negative |
|---|---|---|---|---|---|
| \(k=1\) | 16 | 65,535 | 136 | 2 | 2 |
| \(k=2\) | 21 | 2,097,151 | 2,272 | 21 | 20 |
| \(k=3\) | 24 | 16,777,215 | 13,470 | 54 | 52 |
| \(k=4\) | 25 | 33,554,431 | 11,553 | 73 | 73 |
| total | | 52,494,332 | 27,431 | 150 | 147 |

The \(S_4\times S_4\) class supplies two further records, both certified in
[`certificates/gate_b_n8_k4_rational_v1.json`](certificates/gate_b_n8_k4_rational_v1.json).

**The lowest negative defect certified directly on eight coordinates is \(4/9\).**
The family \(\mathcal T\) of fifteen
rows \(\{0,1,2,60,61,62,64,67,124,128,131,188,193,194,255\}\) has all eight
degrees at the cap \(6\), incidence \(48\ge R_{15}=30\), exactly \(100\) of its
\(225\) ordered joins missing, and

\[
\varepsilon_\vee(\mathcal T)=\frac49=0.4444\ldots,\qquad
A_+(\mathcal T)\le-\frac1{2100}.
\]

This improves the preceding direct record \(139/245\); Proposition 26 later
uses sharp order reweighting to reach \(14/45\). Two features of the direct
family are worth naming. It is *tiny*: low defect turned out not to need a
large family, and
by Lemma 6 a family with defect \(\varepsilon\) needs only
\(m\ge\sqrt{2/\varepsilon}=2.12\) at \(\varepsilon=4/9\), so nothing forced size
here. And it is exactly the shape Corollary 12 demands: four of its eight
coordinates coincide, and it contains \([8]\) itself, well above the dominant-set
threshold \(\frac85\bar s=128/25=5.12\).

**The growth constant improves without leaving the separating class.** The
separating 75-row family \(\mathcal S\) with all degrees at the cap \(30\) and
incidence \(240\ge R_{75}=234\) has

\[
\varepsilon_\vee(\mathcal S)=\frac{284}{375},\qquad
A_+(\mathcal S)\le-\frac{177}{5000},
\]

so its powers give slope \(177/40000=0.004425\) against the published
\(1/250=0.004\): an improvement of \(10.6\%\) using a fully normalized family, so
the growth improvement of Proposition 16 does not depend on admitting
non-separating bases. Its certified ratio \(531/11360=0.0467\ldots\) is the best
among separating families. The cloned base still holds the overall slope record
at \(3/640=0.0046875\).

**The screen is exact at \(k=1\) and \(k=4\) and wrong at \(k=2\) and \(k=3\).**
That is consistent with the diagnosis in Proposition 15: block-pattern orders
represent the objective exactly when the family's automorphism group is the block
group. Of the 150 screen-reported negatives, 147 are certified, four are exactly
non-negative, and one — reported non-negative — is exactly negative.

### Proposition 18: cloning is defect-free; its objective direction is coordinate-dependent

Let \(\mathcal F'\) be \(\mathcal F\) with coordinate \(i\) duplicated. Each row
of \(\mathcal F'\) is a row of \(\mathcal F\) with one bit repeated, so the map
is a bijection commuting with union. Hence, **exactly**:

* \(m\) is unchanged;
* \(\varepsilon_\vee\) is unchanged, since \(X\vee Y\) leaves the family before
  cloning exactly when it leaves after;
* every old degree is unchanged and the clone's degree equals its twin's, so the
  cap still holds;
* the incidence rises by \(\deg_i\) while \(R_m\) depends only on \(m\).

Thus cloning preserves admissibility unconditionally and cannot move the
defect. **Its effect on \(A_+\), however, is not sign-definite.** The earlier
claim that every clone lowers \(Q\), \(C_+\), and \(A_+\) was false. On the
five-coordinate core of \(\mathcal T=\)`n8tiny`, one clone of coordinate 2
lowers the exact upper endpoint from \(+0.010695694102\) to
\(+0.003478204249\), but one clone of any coordinate in
\(\{0,1,3,4\}\) raises it to \(+0.012500066566\).

For the useful coordinate 2, the finite ladder is
([`audit_clone_saturation.py`](audit_clone_saturation.py), verdict
`CHOSEN_CLONE_IS_DEFECT_FREE_AND_FINITE_LADDER_SHRINKS`):

| added clones | \(n\) | \(A_+\le\) | \(\Delta\) | \(\varepsilon_\vee\) |
|---|---|---|---|---|
| 0 | 5 | \(+0.010695694\) | | \(4/9\) |
| 1 | 6 | \(+0.003478204\) | \(-0.007217490\) | \(4/9\) |
| 2 | 7 | \(+0.000772362\) | \(-0.002705842\) | \(4/9\) |
| 3 | 8 | \(-0.000478465\) | \(-0.001250827\) | \(4/9\) |
| 4 | 9 | \(-0.001137229\) | \(-0.000658764\) | \(4/9\) |

This is how the record low-defect negative witness exists: collapsing its four
identical coordinates leaves an admissible five-coordinate family with the same
size 15 and defect \(4/9\) but positive objective. By contrast,
`n8clone_lo` and `n8clone_hi` collapse to inadmissible seven-coordinate
families (incidence 196 below \(R_{70}=215\)), so those clones buy feasibility.

The finite ratios \(0.3749,0.4623,0.5267\) do **not** prove geometric
convergence. Proposition 21 below gives the exact law: with total clone
multiplicity \(r\), the probability that another coordinate precedes the first
clone is \(O(1/r)\). For coordinate 2 the exact limit is
\(A_+\le-0.002338401508\); the total gain is bounded, but convergence is
algebraic rather than geometric.

### Lemma 19: every defect-free coordinate extension is join-consistent

Let \(\mathcal U\subseteq\mathcal F\), add one new coordinate \(*\), and map
\[
\phi_{\mathcal U}(A)=
\begin{cases}
A\cup\{*\},&A\in\mathcal U,\\
A,&A\notin\mathcal U.
\end{cases}
\]
This bijection never repairs a failed union: if \(A\cup B\notin\mathcal F\),
the projection of \(\phi(A)\cup\phi(B)\) onto the old ground set is still
\(A\cup B\), whereas every image row projects into \(\mathcal F\). For a
successful union \(Z=A\cup B\in\mathcal F\), the extended union belongs to the
image exactly when
\[
\mathbf 1_{\mathcal U}(Z)
=\mathbf 1_{\mathcal U}(A)\vee\mathbf 1_{\mathcal U}(B).
\tag{19}
\]
Consequently the defect never decreases, and it is preserved exactly when (19)
holds for every successful ordered pair. Call such a row subset
**join-consistent**.

Equivalently, \(\mathcal U\) is an up-set inside
\((\mathcal F,\subseteq)\), and every successful join landing in
\(\mathcal U\) has at least one part in \(\mathcal U\). Indeed, setting one part
below the other in (19) gives the up-set condition; that condition proves the
forward implication in (19), and join-primality proves the reverse one.

If the old family is admissible and
\(|\mathcal U|\le\lfloor2m/5\rfloor\), the extension is admissible: all old
degrees are fixed, the new degree is \(|\mathcal U|\), and incidence increases
against the unchanged \(R_m\). Thus (19) classifies every defect-free admissible
one-coordinate extension.

There is a useful distinction. A global Boolean function \(g\) commutes with
union on every family iff \(g\) is constant or an OR of input coordinates:
when \(g(0)=0\),
\(g(x)=\bigvee_{i:x_i=1}g(e_i)\), while \(g(0)=1\) forces constant one.
Constant zero is inert and constant one violates the cap. But (19) constrains
only successful joins of the particular family, so defective families can
admit more extensions than global OR-columns.

### Proposition 20: the record core has one useful extension direction

For the five-coordinate core
\[
\mathcal K=(0,1,2,4,5,6,8,11,12,16,19,20,25,26,31)
\]
of `n8tiny`, exhaustive evaluation of all \(2^{15}=32768\) row subsets finds
exactly 57 join-consistent subsets. The cap \(6\) leaves exactly six nonzero
usable columns: the five original coordinate columns, each of weight 6, and
\[
\mathcal U^\star=\{11,19,25,26,31\}
=\{A\in\mathcal K:|A|\ge3\},
\]
of weight 5. The latter is not induced by any global join-homomorphism.
`classify_defect_free_extensions.py` checks the equivalence in Lemma 19 with
zero mismatches and evaluates every usable extension exactly:

* \(\mathcal K\) has \(A_+\le+0.010695694102\);
* cloning coordinate 2 lowers this to \(+0.003478204249\);
* cloning any of \(0,1,3,4\) raises it to \(+0.012500066566\);
* adding \(\mathcal U^\star\) raises it to \(+0.114769081463\).

Thus coordinate 2 is the **only** improving defect-free one-step extension of
this core. Concentrating three clones there gives the certified negative
\(-0.000478465\) witness at unchanged defect \(4/9\), while spreading one clone
over each of coordinates \(0,1,2\) leaves the exact objective positive,
\(+0.006494611\). This is an exact classification for \(\mathcal K\), not a
claim that nonclone extensions are useless for every family.

### Proposition 21: cloning only reweights fixed orders

Replace each original coordinate \(i\) by \(r_i\ge1\) identical copies. In an
order of the cloned coordinates, keep the first copy in every clone class and
delete the others. Once the first copy has been exposed, every later copy is
deterministic in both rows, contributes zero entropy, has one transition of
mass one, and leaves the row fibers unchanged. Therefore, fixed order by fixed
order,
\[
Q_{\widetilde\pi}(\widetilde{\mathcal F})=Q_{\pi}(\mathcal F),
\qquad
C_{+,\widetilde\pi}(\widetilde{\mathcal F})=C_{+,\pi}(\mathcal F),
\qquad
A_{+,\widetilde\pi}(\widetilde{\mathcal F})=A_{+,\pi}(\mathcal F),
\]
where \(\pi\) is the collapsed original order. Averaging merely changes the
probability distribution on the \(n!\) original values, so every cloned
representation lies in their convex hull:
\[
\min_\pi A_{+,\pi}(\mathcal F)
\ \le\ A_+(\widetilde{\mathcal F})\
\le\ \max_\pi A_{+,\pi}(\mathcal F).
\tag{20}
\]

For one coordinate with total multiplicity \(r\), let \(M_j\) average the
original fixed-order objective over orders in which it has rank \(j\).
The first clone has exact rank law
\[
\Pr_r[j]
=\frac{\binom{n+r-j-2}{r-1}}{\binom{n+r-1}{r}},
\qquad
A_r=\sum_{j=0}^{n-1}\Pr_r[j]M_j,
\qquad A_r\longrightarrow M_0.
\tag{21}
\]
In particular \(\Pr_r[j>0]=(n-1)/(n+r-1)\), so the generic convergence scale is
\(O(1/r)\), not geometric.

This exactly explains `n8tiny`: its minimum fixed-order lower endpoint is
\(-0.003215773322\), coordinate 2 has limit
\(-0.002338401508\), and (21) crosses zero at total multiplicity four.

It also closes cloning on the new normalized low-defect core
\[
\mathcal R=(0,2,4,6,8,9,16,22,31,32,96,105,112,125,127)
\subseteq2^{[7]}.
\]
Here \(m=15\), degrees are \((5,5,6,6,6,6,5)\), incidence is
\(39\ge R_{15}=30\), \([7]\in\mathcal R\), and
\(\varepsilon_\vee(\mathcal R)=14/45<2/5\). Its exact order average is
\[
A_+(\mathcal R)\in
[0.305904110097,\ 0.305904110098],
\]
and **every one of its 5040 fixed-order lower endpoints is at least
\(0.145914214376\)**. Hence (20) proves that no clone multiplicities, on one
coordinate or many, can make it negative.

Lemma 19 finds 74 join-consistent subsets of \(\mathcal R\), with eight usable
under the cap: the seven old columns and exactly one new column. Adding that
column preserves defect and admissibility, gives
\[
A_+\in[0.295857317023,\ 0.295857317024],
\]
and every one of its 40320 fixed-order lower endpoints is at least
\(0.122963869853\). A defect-free extension preserves the successful-join
relation on row identities, so its join-consistent subsets are the same eight;
after the one new column has been used, every further extension is a clone.
The two convex-hull certificates therefore prove that **every finite sequence of
defect-free extensions of \(\mathcal R\) remains positive**. This is
`audit_clone_limit.py`, verdict
`DEFECT_FREE_EXTENSIONS_CLASSIFIED_AND_SUB_TWO_FIFTHS_CORE_CANNOT_BE_RESCUED`.

### Proposition 22: no cap-only dominant-set floor exists; Reimer is essential

Chase--Lovett [REPORTED, arXiv:2211.11689v1, Example 1.4] construct
\[
\mathcal F_n=
\{A:|A|=\psi n+n^{2/3}\}
\cup\{A:|A|\ge(1-\psi)n\},
\qquad \psi=\frac{3-\sqrt5}{2}.
\]
They prove \(\varepsilon_\vee(\mathcal F_n)=o(1)\) and every frequency is at
most \(\psi+o(1)<2/5\). The construction is active and separating for large
\(n\), and it contains \([n]\). Therefore a positive defect floor does **not**
follow from the cap, under either the displayed cap-only convention or the
normalized search convention, even in the presence of a dominant set.

The construction is not Gate B admissible. Its first layer dominates, so
\[
\frac1n\log_2|\mathcal F_n|\longrightarrow h(\psi),
\qquad
\frac{\bar s}{n}\longrightarrow\psi,
\]
whereas Reimer demands
\(\bar s/n\ge h(\psi)/2\). Exact rational entropy enclosures give
\[
\frac{h(\psi)}2-\psi>0.0977433528.
\]
Thus it misses Reimer by a linear amount. Any defect floor for the actual Gate
B class must exploit Reimer essentially; frequency, pair defect, normalization,
and a dominant set are jointly insufficient.

The same audit locates the entropy-method obstruction finitely. A subset-DP
computes \(\max_\pi Q_\pi\) without enumerating \(n!\) orders. Six of the ten
registered negative bases satisfy
\(\max_\pi Q_\pi<\log_2m\), so the iid Gilmer chain fails at every order; the
widest exact margin is at least \(0.0393134578\). All ten registered bases have
maximum frequency exactly \(2/5\). Cambie's reported two-atom obstruction is
enclosed at
\[
c^\star\in
[0.382345533366702721,\ 0.382345533366702722],
\]
but its matching lower verification remains **OPEN / COMPUTATIONAL-EVIDENCE**
in the broader `math/uc` audit and is not promoted here.

The computational attack was then run in the order the proof dictates. A
deterministic defect-only search at \(n=9,10\), covering 13 and 15 sizes respectively including maxima \(145,255\),
found admissible families below \(2/5\) at
\[
\frac{14}{45}\quad(n=9,m=15),\qquad
\frac{8}{25}\quad(n=10,m=15),
\]
so the former \(692/2025\) searched floor was not structural. The first value
collapses to the normalized core \(\mathcal R\) treated in Proposition 21 and is
exactly positive under every defect-free extension.

A second search held \(\varepsilon_\vee<2/5\) and \([n]\in\mathcal F\) as hard
constraints and minimized the full-order iid cost \(Q\). Only one finalist
passed the necessary test
\((1-\alpha)Q-\log_2m\le0\): an \(n=9,m=20\) family at defect \(19/50\).
Its Q-only lower endpoint is \(-0.003072867378\), but the sampled Bellman term
raises the sampled objective to \(+0.168611231079\): a \(+0.171684\ldots\)
contribution, with sampled \(C_+\approx4.82\). Thus \(C_+\), not \(Q\), is the
binding term on this survivor.
After deleting its one duplicate column, the normalized eight-coordinate core
has exact
\[
A_+\in[0.172732849547,\ 0.172732849548],
\]
and all 40320 fixed-order lower endpoints are at least \(0.042394032883\).
Proposition 21 therefore excludes every cloning pattern on that survivor.

Finally, a sampled-order Bellman walk optimized the full \(A_+\) under the same
hard defect and dominant-set constraints. Its best discovery value remained
\(+0.128833928729\) at defect \(86/225\), more than twelve times the
predeclared \(+0.01\) threshold for exact clone follow-up. This last value is
**DISCOVERY ONLY**; the exact conclusions are the defect, the full-order
\(Q\) enclosures, and the two all-order convex-hull audits above. No certified
negative family below \(2/5\) was found.

The search implication is methodological: use full-order \(Q\) as a necessary
rejection gate, but rank future walks by a cheap one-sided-cost proxy or sampled
full \(A_+\); a \(Q\)-only ranking preferentially retains candidates that the
Bellman audit then kills.


### Proposition 23: Reimer and the dominant set are automatic at zero defect

Let \(\mathcal F\) satisfy the \(2/5\) cap and
\(\varepsilon_\vee(\mathcal F)=0\). Then \(\mathcal F\) is union closed, so
Reimer's average-set-size theorem [CITED-DEPENDENCY: D. Reimer, *An Average Set
Size Theorem*, Combinatorics, Probability and Computing 12(1) (2003), 89--93,
doi:10.1017/S0963548302005230] gives
\[
\bar s\ge\frac12\log_2|\mathcal F|.
\]
Thus the repository's Reimer admissibility condition adds no restriction at
zero defect.

The dominant-set condition is automatic there as well. The finite union
\(T=\bigcup_{A\in\mathcal F}A\) belongs to \(\mathcal F\). On the active
support, \(T=[n]\), while the cap gives
\[
\bar s=\sum_{i=1}^n\frac{\deg_i}{m}\le\frac{2n}{5}.
\]
Consequently
\[
|T|=n\ge\frac52\bar s\ge\frac85\bar s.
\]
The second inequality is strict when the active support is nonempty; equality
handles \(\mathcal F=\{\varnothing\}\). Thus \(T\) meets the dominant-set
threshold required by Corollary 12 in every case.

This remains true under the normalized search convention. Delete every zero
column and all but one column in each duplicate class. Projection commutes with
union and is injective on the row family: a deleted zero column carries no
information, and every deleted duplicate is recovered from its retained
representative. The projected family has the same \(m\), defect zero, and cap,
is active and separating, and satisfies Reimer again by the cited theorem.

Hence any positive defect floor obtained from cap, Reimer, a dominant set, and
optionally normalization would prove that every nontrivial union-closed family
has a frequency strictly above \(2/5\). Reimer is essential for positive defect,
by Proposition 22, but becomes automatic exactly at the endpoint. This is the
precise reason a floor reaching \(\varepsilon_\vee=0\) is the \(2/5\) frequency
problem rather than a perturbative consequence of Reimer.

### Proposition 24: Reimer has no pair-defect-continuous extension

Put
\[
\Delta_R(\mathcal F)
:=\frac12\log_2|\mathcal F|-\bar s.
\]
There is no function \(g(t)\to0\) as \(t\to0^+\) for which every active,
separating cap-\(2/5\) family containing its full set satisfies
\[
\Delta_R(\mathcal F)\le n\,g(\varepsilon_\vee(\mathcal F)).
\tag{22}
\]

Indeed, the Chase--Lovett sequence in Proposition 22 has all those structural
properties and \(\varepsilon_\vee(\mathcal F_n)=o(1)\), but
\[
\frac{\Delta_R(\mathcal F_n)}n
\longrightarrow \frac{h(\psi)}2-\psi>0.0977433528.
\]
The convergence statements are [REPORTED] from arXiv:2211.11689v1, Example
1.4; the strict decimal lower bound is the exact rational entropy enclosure in
`barrier_sawin_cambie.py`. Equation (22) would make the left side tend to zero,
a contradiction.

Thus approximate union closure in the ordered-pair sense cannot recover even an
approximate Reimer theorem with vanishing normalized error. This does not rule
out a bound on families that already satisfy Reimer; it rules out manufacturing
that hypothesis continuously from \(\varepsilon_\vee\), even after adding the
cap, normalization, and a full set.

### Proposition 25: cloning reaches exactly the fixed-order extrema

Replace original coordinate \(i\) by \(r_i\ge1\) labeled identical copies.
Proposition 21 proves that deleting every copy after the first occurrence in a
cloned order leaves the fixed-order objective unchanged. The distribution of
the retained original order is also exact. If
\(\pi=(\pi_1,\ldots,\pi_n)\), then
\[
\Pr_{\mathbf r}(\pi)
=\prod_{k=1}^n
\frac{r_{\pi_k}}{\sum_{j=k}^n r_{\pi_j}}.
\tag{23}
\]
To see this, give every labeled copy an independent rate-one exponential clock.
Their clock order is uniform, the first clock in class \(i\) is exponential of
rate \(r_i\), and the exponential race plus memorylessness gives each factor in
(23). Therefore
\[
A_+(\widetilde{\mathcal F}_{\mathbf r})
=\sum_{\pi\in S_n}\Pr_{\mathbf r}(\pi)A_{+,\pi}(\mathcal F).
\tag{24}
\]

Conversely, fix a target order \(\pi\) and an integer \(R\ge2\), and set
\[
r_{\pi_k}=R^{\,n-k}.
\]
At each of the first \(n-1\) stages the chosen weight is at least a
\((1-1/R)\)-fraction of the remaining weight. Hence
\[
\Pr_{\mathbf r}(\pi)
\ge(1-1/R)^{n-1}
\ge1-\frac{n-1}{R}
\longrightarrow1.
\tag{25}
\]
Combining (20), (24), and (25) gives the sharp closure
\[
\inf_{\mathbf r\in\mathbb Z_{\ge1}^n}
A_+(\widetilde{\mathcal F}_{\mathbf r})
=\min_{\pi\in S_n}A_{+,\pi}(\mathcal F),
\qquad
\sup_{\mathbf r\in\mathbb Z_{\ge1}^n}
A_+(\widetilde{\mathcal F}_{\mathbf r})
=\max_{\pi\in S_n}A_{+,\pi}(\mathcal F).
\tag{26}
\]
These are infimum and supremum statements; finite attainment is not asserted.
Sign reachability is finite and exact. If
\(a_{\min}:=\min_\pi A_{+,\pi}<0\) and
\(a_{\max}:=\max_\pi A_{+,\pi}\), then (24)--(25) give
\[
A_+(\widetilde{\mathcal F}_{\mathbf r})
\le a_{\min}+\frac{n-1}{R}(a_{\max}-a_{\min})<0
\]
whenever
\[
R>\frac{(n-1)(a_{\max}-a_{\min})}{-a_{\min}}.
\]
If \(a_{\min}\ge0\), every convex combination is nonnegative. Therefore a
family can be made negative by a finite clone multiset **if and only if** one
of its original fixed orders is negative.

`audit_clone_limit.py` implements (23)--(26) in exact rational arithmetic;
`test_gate_b.py` checks (23) against uniform multiset-order enumeration. The
audit constructs a finite multicoordinate witness whenever a negative order is
certified. For the earlier defect-minimizing \(14/45\) core \(\mathcal R\), all
5040 fixed orders remain positive, so the same criterion is a complete
impossibility theorem, not merely an asymptotic budget estimate.

Cloning preserves \(m\), cap/Reimer admissibility, and
\(\varepsilon_\vee\), but a nontrivial clone creates duplicate columns.
Accordingly the constructive half applies to the displayed cap/Reimer
supremum, while the normalized convention records the result as
non-separating. The impossibility half applies to every cloning attempt under
either reporting convention.


### Proposition 26: the certified negative defect reaches \(14/45<2/5\)

The minimum-fixed-order row search found the normalized core
\[
\mathcal V=(0,1,2,4,5,8,10,43,64,190,192,193,245,254,255)
\subseteq2^{[8]}.
\]
It has \(m=15\), degrees \((6,6,6,6,4,5,6,6)\), incidence
\(45\ge R_{15}=30\), contains \([8]\), and has exact defect
\[
\varepsilon_\vee(\mathcal V)=\frac{14}{45}.
\]
Its uniform order average is positive,
\[
A_+(\mathcal V)\in[0.230029313333,\ 0.230029313334],
\]
but the fixed order
\(\pi=(6,1,2,0,3,4,5,7)\) has the exact rational enclosure
\[
A_{+,\pi}(\mathcal V)\in
\left[
-\frac{280488497585020501936812613}
       {79228162514264337593543950336},
-\frac{140244248792510250968406301}
       {39614081257132168796771975168}
\right]
\subset(-\infty,0).
\tag{27}
\]
The order was selected by a sampled search, but (27) is an exact two-sided
enclosure from `verify_gate_b_rational.py`; the sampled value is not used after
selection.

Apply Proposition 25 with ratio \(R=41\) and coordinate multiplicities
\[
\mathbf r=(2825761,\ 4750104241,\ 115856201,\ 68921,\
1681,\ 41,\ 194754273881,\ 1).
\]
The resulting symbolic cloned family \(\widetilde{\mathcal V}\) has dimension
\(199623130728\), the same \(15\) rows and defect \(14/45\), maximum degree
\(6=\lfloor2m/5\rfloor\), and incidence
\(1197738780965\ge R_{15}\). Exact reweighting of all \(8!=40320\) fixed-order
enclosures gives
\[
A_+(\widetilde{\mathcal V})\in
\left[
-\frac{3339001091075785382516669}
       {39614081257132168796771975168},
-\frac{3339001091075785382516663}
       {39614081257132168796771975168}
\right]
\subset(-\infty,0).
\tag{28}
\]
Thus the lowest certified negative defect in the displayed cap/Reimer class
drops from \(4/9\) to
\[
\boxed{\frac{14}{45}=0.3111\ldots<\frac25}.
\]
The core \(\mathcal V\) is active and separating, but
\(\widetilde{\mathcal V}\) is not separating because it has repeated columns.
Therefore the separating/normalized negative frontier remains \(1144/1875\);
the new \(14/45\) record is explicitly under the displayed convention.

The authoritative artifact is
`certificates/gate_b_subtwofifths_clone_rational_v1.json`, verdict
`PROVED_SUB_TWO_FIFTHS_NEGATIVE_EXACT_RATIONAL`. Its frozen targets are the
exact rationals \(A_+<0\) and \(\varepsilon_\vee<2/5\), not float values. The
row search is heuristic coverage; every property consumed in (27)--(28) is
rechecked exactly.

## Scope and stronger interpretation

This theorem resolves the repository's displayed universal supremum. The
constructed defects tend to one, not zero. Therefore it does **not** resolve a
different local-stability question restricted to sequences with
\(\varepsilon_\vee\to0\). The prior plan's phrase "near-UC adversarial search"
was a method proposal; no such restriction appeared in the Gate B definition.
The distinction is recorded explicitly rather than silently changing Gate B.

For that local question, Proposition 26 crosses the formerly hard \(2/5\)
barrier, but it does not determine \(c_{\rm loc}\): its defect is the fixed
positive value \(14/45\), not a sequence tending to zero. Proposition 25 now
gives the exact mechanism test. A further reduction must find row-changing
cores with still smaller defect and a negative fixed order, or prove a
Reimer-essential floor. At the endpoint, Proposition 23 shows that a positive
floor would yield the \(2/5\) frequency theorem; the present work neither
assumes nor manufactures that theorem.

The result also does not produce a union-closed counterexample or prove
Frankl's conjecture. It proves that every fixed scalar multiple of the ordered
pair closure defect is structurally incapable of repairing this particular
\(A_+\) relaxation over all cap/Reimer families.
