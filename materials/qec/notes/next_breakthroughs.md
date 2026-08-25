# Where the next substantially stronger results are — ranked, with first experiments

Written 2026-08-21 after Theorems J-G/J-H/J-I/J-J (EXP-053/054). Ranking is by
*expected information gain × payoff / cost*, not by appeal. Every entry names a
falsifiable prediction and the cheapest decisive first experiment, so a wrong
guess dies fast.

State that makes these live now:

* the single-row X-collapse channel is **completely classified** by one ideal
  ($I=\operatorname{Ann}_R(A,B)$): nilpotent → every direction collapses,
  idempotent → none, otherwise mixed (Theorem J-G);
* the collapse **cannot happen at all on odd$\times$odd lattices** (J-G1);
* it **cannot be partial** for weight-$\le3$ polynomials (J-I) — the entire
  published family — but is partial at weight 4 (J-J);
* on the Z-side the trade law $k_Q=k_P-T$ is exact on 138/139 certified parents,
  and 63 parents (incl. the Gross code) are family-closed.

---

## A. Odd-lattice PBB search — exact closure through $n=234$ (EXP-055--067)

**Status 2026-08-21: the "unexplored region" premise was FALSE and is retracted.**
J-G1 does say the collapse hazard is structurally absent when $\ell,m$ are both
odd, and our own 202-parent catalogue really does live at $m\in\{3,6\}$. But the
claim that "$[[90,8,10]]$ on $(15,3)$ is the only odd$\times$odd BB instance in
the literature" was wrong: the validation battery now contains **29 sourced
odd$\times$odd BB instances** in Wang–Mueller (arXiv:2408.10001v4),
Eberhardt–Steffan (arXiv:2407.03973v1), Postema–Kokkelmans
(arXiv:2502.17052v4), Bravyi et al. (arXiv:2308.07915), and Liang et al.
(arXiv:2503.03827v3). The later audit also missed Liang et al.'s two exact
$[[234,8,18]]$ Table III rows; that omission is corrected. The region is
*searched*; only our own catalogue avoided it.

**What was also already published:** the odd-lattice rate law itself. $R$
semisimple $\Rightarrow k=2\cdot\#\{\text{common roots}\}$ is
Panteleev–Kalachev (arXiv:1904.02703) Prop. 1 in the cyclic case,
Lin–Pryadko (arXiv:2306.16400) Eq. (47), Postema–Kokkelmans Thm. 2.6, and
Eberhardt–Steffan Cor. 2.11–2.12 ("if $\ell$ and $m$ are odd, all BB codes are
principal"). Claim it as reproduction, never as novelty.

**What survived and is now done (EXP-055--067).**

1. *A certified solver-free distance ceiling* (the exact pole isomorphism is
   published; its explicit use as a general-BB rejection oracle was not found).
   If $I=\operatorname{Ann}_{\rm left}(a,b)$ and $J=\bar I$ is the physical
   right-kernel pole, then $J^2\cong\ker H_X/S_Z$ exactly and
   $d\le\min\{\operatorname{wt}(u):0\ne u\in J\}$. Reduced pole-coset witnesses
   are self-certifying and much tighter. The isomorphism passes 29/29
   literature audits; the ceiling has zero violations on eight locally exact
   rows (all 27 reproduced reported values also pass as a non-certifying
   sanity check).
2. *An exhaustive census*: all $65$ odd lattices with $\ell m\le180$,
   $4.23\times10^9$ weight-$\le3$ pairs, exact $k$ by two independent routes
   (zero mismatches), immunity re-verified ($273$ idempotence tests, zero
   violations).
3. *An exact orbit-class solver*: EXP-056 reduces the $(3,27)$
   $[[162,8,\cdot]]$ target's 255 nonzero classes to 20 verified automorphism
   orbits. Sparse-generator CNFs are UNSAT 20/20 and replay 20/20; the
   even-kernel law excludes 13 and an explicit weight-14 word closes
   $d_X=d_Z=14$. The old 2,100 s CP-SAT minimization timeout is superseded.
4. *An adaptive exact frontier ratchet*: odd-column parity skips odd caps,
   verified BB duality halves sector work, and raw kissat stops at first
   SAT/prior UNSAT. It exactifies $[[170,16,10]]$, seven
   $[[186,10,14]]$, thirteen $[[210,18,8]]$, $[[210,24,4]]$,
   $[[210,14,12]]$, and two $[[210,10,16]]$ classes. A separate $k=8$
   row is dominated by independent weight-16 witnesses; its distance is not
   claimed exact.
5. *A cross-lattice CRT theorem*: all coprime factorizations of fixed
   $N=\ell m$ are explicitly coordinate-permutation equivalent. Thus the
   $(15,7)$, $(21,5)$ and $(35,3)$ screens are one exact problem.
   Monotone reference rebinding then preserves old proofs without re-solving.

**Screen verdict.** Complete through $n=210$: 20 nonempty-frontier lattices,
4,020 symmetry classes / 122,833 represented pairs; all 3,816 classes with an
independently exact reference are dominated (3,364 reduced-pole witnesses,
395 bounded/adaptive CDCL witnesses, 57 exact CP-SAT fallbacks), 204 high-$k$
no-reference, zero survivors or undecided. All 3,759 explicit witnesses are
physically rechecked.

**Method falsifiers.** The exact affine-trellis preflight is dead on the N=105
hard core: best widths 48/53/80/99 exceed the five-million-state gate. The
$H=7$ quotient projection is exact but yields lower bound only 4. Pole
stabilizer $H=5$ characterizes all thirteen easy $[[210,18,8]]$ classes, but
``hard implies $H=1$'' is falsified by hard classes with $H=5$ and $H=7$.

**Exact $n=234$ closure (EXP-066/067).** The two frontiers are exhaustively
represented: $(13,9)$ closes 84/84 classes; on $(39,3)$ the full 576-element
automorphism action compresses 182 hard classes to 30 bundles and transported
witnesses close 28. EXP-067 uniquely maps Liang et al.'s two published
$[[234,8,18]]$ twisted-torus rows to the two residual representatives. For
each, original and block-swapped rooted connected-cluster searches exhaust
cap 16 and replay; even kernel weight and physical weight-18 witnesses give
$d_X=d_Z=18$. Promoting that reference closes all 4,658/4,658 referenced
classes across 22 lattices / 4,862 classes / 150,581 pairs, with 204
no-reference and zero survivor/undecided. The multithreaded solver coordinator
was rejected for FR-033's completion race; only the source-bound single-thread
entrypoint is certifying.

**Next A-step.** Screen the $n=270$ frontiers. First collapse the isomorphic
noncyclic presentations $(15,9)$ and $(45,3)$; then screen their 5,024 classes
once and the cyclic $(27,5)$ lattice's 688 classes. Liang Table III's two
$[[270,8,20]]$ rows are source leads, not local exact references until their
constructors and lower bounds are bound by the current protocol.

## B. Multi-row light channels: replace CP-SAT certificates with a theorem

**Gap.** J-E/J-G decide the *single-row* channel exactly and instantly. The
multi-row channel is decided today by a one-shot CP-SAT query per parent
(EXP-051); it is why `15_6_0256` and `30_6_0289` remain X-undecided, each
carrying a verified weight-8 channel.

**Target theorem.** Generalize J-E from rank-1 dressings $\lambda[C\,D]$ to
rank-$r$: with $\Delta$ the dressing submodule, characterize
$\exists w<d_X(P)$ in $\operatorname{rowspace}(H_X)\setminus\{\text{original
rows}\}$ that becomes logical, in terms of the ideals generated by
$r\times r$ minors of $[C\,D]$ over $R$ — a Nakayama/Fitting-ideal statement.
J-H's local decomposition $R_\chi=F_\chi[G_2]$ is the right coordinate system:
per factor the question becomes a small linear-algebra condition over $F_\chi$.

**Prediction.** The multi-row channel is empty whenever the *Fitting ideal* of
the dressing module is nilpotent — the exact analogue of J-G, testable
immediately against EXP-051's 8 INFEASIBLE certificates and the 2 witnesses.

**First experiment.** Compute the Fitting-ideal invariants for all 10 immune
parents and correlate with EXP-051's verdicts (8 empty / 2 nonempty). Cost:
hours. A clean correlation is strong evidence; a mismatch kills the conjecture
cheaply.

## C. Weight-4 mixed map, and the first *usable* PBB safety subspace

**Why it matters practically.** In a mixed parent the fixed set $S$ is a
*proper nonzero* subspace: perturbation directions inside $S$ provably cannot
demote, those outside provably do. That is the first structure that would let a
designer *use* the PBB freedom safely — pick $[C\,D]$ inside $S$ (dimension
$2\dim I^\infty$, explicitly computable). Today the only mixed instance is a
$[[12,8]]$ toy.

**Target.** The sharp weight/lattice map at weight 4: for which $(\ell,m)$ and
coset patterns does mixed occur, and is there a mixed parent with $k_P>0$ and
certified $d\ge6$ at $n\le180$?

**First experiment.** Weight-4 census with EXP-054's group-dedup algorithm
(per-polynomial invariants, then group$\times$group pair tests) on the catalogue
lattices; $\binom{\ell m-1}{3}$ normalised supports is $\le10^6$ per lattice, so
$(6,6),(9,6),(12,6),(6,9),(4,9)$ are directly feasible and $(15,6),(30,6),
(12,12),(15,12)$ need the cached invariants plus a coset-pattern prefilter
derived from J-H (only patterns with $\ge2$ occupied cosets and $\ge2$ points per
occupied coset can be mixed — a large pruning).

**Prediction from J-H.** Mixed at weight 4 requires a lattice with a nontrivial
$G_2$ *and* a nontrivial $G_{\rm odd}$ (so $(2^s,2^t)$-only and odd-only
lattices are excluded), and each of $A,B$ must be a product
$(1+h)\cdot(\text{odd-part polynomial with a character zero})$ up to
translation. That is a concrete constructive recipe — it can be tested by
building the recipe's outputs and checking `case == mixed` directly.

## D. Rate–distance obstruction from the ideal decomposition

**Status of the naive form: dead, recorded.** For immune parents all logical
classes live in the idempotent ideal $e_IR$ (a 2D cyclic code of GF(2)-dimension
$k_P/2$), which gives $d_X(P)\le\min\operatorname{wt}(e_IR\setminus 0)$. Measured
on all 10 immune parents: the ideal minimum weight is $18$ on $(9,6)$ and $30$ on
$(15,6),(30,6)$ against certified $d_P\in\{6,8,10\}$ — the bound is loose by
$2$–$3\times$ because $\ker H_X$ is strictly larger than $e_IR^2$ (syzygy pairs
$(u,v)$ with $\bar au=\bar bv\neq0$ carry the true minimum). Recorded in
`failed_routes.md`.

**Live version.** Bound the *coset* minimum in $\ker H_X/S_Z$ using the local
factor decomposition: on each factor the code is a small cyclic-code problem
over $F_\chi[G_2]$, and distances of 2D cyclic codes have classical bounds
(BCH-type over the character set). A per-factor bound would give the first
*algebraic* $d$-upper-bound for BB codes that does not require a solver — the
tool the envelope work has been missing at $n=360$.

**First experiment.** For the 10 immune parents and the Gross code, compute the
per-factor character supports and compare classical cyclic bounds against the
certified distances; look for the tightest factor.

## E. The decisive end-to-end item: a fault-tolerant mixed-stabilizer circuit

**Where it stands.** Theorem C6 proves *every* generating set of a catalogue PBB
$[[144,12,12]]$ has a symplectic-weight-$\ge8$ generator, so any one-ancilla
schedule needs $\ge8$ two-qubit layers against the Gross code's 7 — and no
fault-tolerant one-ancilla circuit for the mixed stabilizers was found in the
searched schedule space (EXP-009/016/031/036/039). **No impossibility is
claimed**, and that is the honest gap in the original assignment.

**Two routes, both concrete.** (i) *Construct*: allow depth $\ge8$ with flag
qubits or cat states, verify noiseless semantics, generate the DEM, and compute
or bound the circuit distance (the EXP-040 machinery already does mechanism
distance on flattened DEMs). (ii) *Prove*: a no-go for one-ancilla schedules of
mixed generators within a stated class — the incidence-counting route is already
refuted (ancilla-ancilla parity forwarding), so the lower bound must come from
the symplectic weight structure, i.e. from C6 plus a hook-error argument.

**Payoff.** Route (i) turns the provisional matched-MC pair (Gross
$6.33\times10^{-3}$ vs PBB $2.21\times10^{-2}$, disjoint CIs) into a real
end-to-end comparison. Route (ii) is a publication-grade no-go on its own.

## F. Rare-event estimator to finish the matched comparison

The matched Monte-Carlo grid is 1/10 cells complete because ordinary MC needs
too many shots at $p\le0.002$. A validated splitting / importance-sampling
estimator (with an overlap-regime calibration against ordinary MC and exact
enumeration on small circuits) both finishes the grid and is itself a
publication-grade deliverable under the assignment's criterion 9. Cost: days;
risk: estimator bias must be demonstrated controlled, not assumed.

---

## Ordering and dependencies

```mermaid
graph TD
  A[A. odd-lattice PBB sweep] --> E[E. FT mixed-stabilizer circuit]
  B[B. multi-row Fitting-ideal theorem] --> A
  C[C. weight-4 mixed map + safety subspace] --> A
  D[D. per-factor distance bounds] --> A
  D --> G[n=360 envelope closure]
  F[F. rare-event estimator] --> E
```

A's first bounded phase is closed through $n=162$; its next extension needs a
new exact reference, not another blind sweep. B and C are independent theory
tracks; D unlocks the remaining $n=360$ envelope rows. E remains the decisive
end-to-end item, and F is what makes E affordable.
