# General-size GRS circuit transport, C2 at the ambient top, C1/C3 layers
## gate H-GRS-GENERAL-SIZE, run `20260902T113614Z_a999d1b0_5daa598b6a04`

Agent `RsPe3dH2`, 2026-09-02 UTC. Preregistration:
`prereg/H_GRS_GENERAL_SIZE_PREREG_2026-09-02.md`, path-scoped source commit
`af2274660a3b0efcfb38c0f4e7499292fee27894`, SHA-256
`94bbc88be7d500b303dfd95656953045a19f7bda9be52c835bb454acd48210b9`. The in-run
`pre_statement.md` is byte-identical to that source (`cmp` and SHA-256 both
checked); the preregistration was not amended. No parameter-dependent
computation happened before the preregistration commit.

## 1. Theorem T (transport) — verified with obligations (a)–(c)

For GRS/Vandermonde factors $a_x=(1,x,\dots,x^{r_A-1})^{\mathsf T}$ and
$b_y=(1,y,\dots,y^{r_B-1})^{\mathsf T}$ at distinct evaluation points, with
arbitrary nonzero column multipliers, the product column at $(x,y)$ is the
row-major Kronecker evaluation $h(x,y)=a_x\otimes b_y=(x^i y^j)$,
$0\le i<r_A$, $0\le j<r_B$, of dimension $r_A r_B$.

In-run verification at GF(7), GF(11), GF(13):

- **(a) Column form / multipliers rank-free.** Rescaling every product column
  by distinct nonzero scalars leaves the rank unchanged
  (`stageA.multipliers_rank_free`).
- **(b) Kronecker-vs-slice identity.** The matrix of
  $\Gamma\mapsto A\Gamma B^{\mathsf T}$ in row-major vec order is exactly the
  Kronecker evaluation matrix: its rows equal $h(x_u,y_j)$, and
  $(A\otimes B)\operatorname{vec}(\Gamma)=0\iff A\Gamma B^{\mathsf T}=0$ on
  structured sparse test systems at all three primes
  (`stageA.obligation_b_identity.p{7,11,13}`).
- **(c) Minimality transport.** All-support circuit tests used actual
  Kronecker columns with two independent rank algorithms (batch RREF and
  incremental basis, cross-checked on every computation through `dual_rank`);
  the deletion predicate (C) was exercised on every counted support. Ties
  collapse as predicted: the diagonal 6-tie at (2,4) has rank 5, not 6
  (`stageA.tie_collapses_rank`), confirming the all-distinct stratum is where
  the size-$(r_A+r_B)$ law lives.

## 2. C1 — forward law and per-prime PGL equalities

**C1 (forward, proved).** For all-distinct $S$ with $|S|=r_A+r_B$ whose pairs
lie on one nondegenerate Möbius graph $y=-\dfrac{c_2x+c_0}{c_3x+c_1}$,
$\Delta=c_0c_3-c_1c_2\ne0$, $S$ is dependent with rank $r_A+r_B-1$. Proof
sketch registered in the prereg: the graph equation is a nonzero bidegree-
$(1,r_B-1)$-and-$(r_A-1,1)$-span relation; in exact arithmetic the run
verified dependence at every anchor (all 15 sweeps below show
`dependent = circuits = pgl_sum`).

**C1 per-prime equality at the anchors (measured, support-set exact).** For
points $\{1,\dots,n\}\bmod p$, all-distinct $k$-sets, $k=r_A+r_B$:

| $(r_A,r_B)$ | $n$ | $k$ | GF(7) | GF(11) | GF(13) |
|---|---:|---:|---:|---:|---:|
| (3,3) | 6 | 6 | 12 | 4 | 2 |
| (3,3) | 7 | 6 | 588 | 72 | 42 |
| (2,4) | 6 | 6 | 12 | 4 | 2 |
| (2,4) | 7 | 6 | 588 | 72 | 42 |
| (3,4) | 7 | 7 | 42 | 2 | 2 |

At every cell the direct exhaustive sweep, the Möbius-graph support set, and
the PGL sum $\sum_M\binom{|D_M|}{k}$ agree as SETS (not merely counts):
`stageC.*.eq_a_support_set` and `stageC.*.eq_b_count` all pass. Sub-$k$
emptiness at the anchors: all-distinct dependent 5-sets at (2,4) are 0 and
6-sets at (3,4) are 0 over all three primes (`stageC.*.k5_empty`,
`stageC.34.*.k6_empty`). PGL enumeration sizes were asserted
$|\mathrm{PGL}(2,p)|=p(p^2-1)$: 336 / 1320 / 2184.

**(2,2) roots.** The (2,4), (3,3), and (3,4) cells above carry the owner
anchors: $(2,4)$ $n=6\to2$, $n=7\to42$ at GF(13); $(3,4,7)\to2$ GF(13)/GF(11)
and $42$ at GF(7) (where the residues of $1..7$ are $1..6,0$ and the size-7
all-distinct layer is nonempty; GF(7) is the one prime where $n=7$ does not
avoid the folded residue).

**SCOPED (not certified):** the converse of C1 (every all-distinct dependent
$(r_A+r_B)$-set lies on a single Möbius graph) — true on every enumerated
anchor instance, unproved in general; general emptiness below $r_A+r_B$
beyond these anchors (slice-device boundary recorded in the prereg). Neither
gates any verdict here.

## 3. C2 — the two top sizes (proved; censuses pinned)

**C2 (i), size $r_Ar_B$.** $S$ is a circuit iff $E_S$ has corank 1 with every
one-column deletion of full rank $r_Ar_B-1$; pure linear algebra, verified
in-run on every counted support by explicit deletion witnesses.

**C2 (ii), size $r_Ar_B+1$.** Dependence is automatic (ambient cap); $S$ is a
circuit iff every $(r_Ar_B)$-subset is independent, equivalently all
$|S|$ deletion minors are nonzero. Re-expresses gate 9's predicate (14)
exactly.

Censuses (exhaustive, per prime), all matching the preregistered pins:

| grid | sizes 9 / 10 circuits | GF(7) | GF(11) | GF(13) |
|---|---|---:|---:|---:|
| $4\times4$ | 9 | 0 | 0 | 0 |
| $4\times4$ | 10 | 16 | 16 | 16 |
| $5\times5$ | 9 | 36378 | — | — |
| $5\times5$ | 10 | 181960 | — | — |

Gate 9's certified $5\times5$ GF(7) counts 36378 / 181960 are reproduced
exactly by the C2 predicates — the re-expression claimed in the prereg.

## 4. C3 — pencil layer at size $r_Ar_B-1$ ((3,3) bootstrap certified)

**C3 (proved direction).** Two independent bidegree-$(2,2)$ forms without a
common component meet in $(2,2)\cdot(2,2)=8$ points with one residual; the
size-8 support is a genuine circuit of $A\otimes B$ at (3,3). Verified by the
gate-10 engine (Sylvester resultants both ways plus content GCDs, with
dual-route agreement asserted inside `common_component_info`).

**Bootstrap census, $n=8$ complete matching sweeps** (all $8!=40320$
bijections per prime):

| prime | rank-7 pencils | complete intersections | common-$(1,1)$ noncircuits |
|---:|---:|---:|---:|
| GF(11) | 2192 | 560 | 1632 |
| GF(13) | 1344 | 416 | 928 |

These match gate 10's certified counts exactly
(`20260902T041030Z_64e0c4ef_b5de16e49f8a`). Every CI support's resultant
polynomials were additionally compared against node polynomials where the
route applied.

**SCOPED (not certified):** the no-fixed-component characterization of the
whole slice — the converse fails at (3,3) (gate-10 shared-component pencils);
the control `controls.shared_component_reject` plants that failure.

## 5. Controls (both directions; final count 103)

Preregistered §7 battery, all passing in the final complete run:

1. **Pure-tensor ACCEPT/REJECT** — verified size-6 (2,4) Möbius circuit
   accepted; 128 corrupted variants tried, first coordinate shift
   (`delta=1`) raises rank and is rejected
   (`controls.tensor_accept`, `controls.corrupted_reject`). (2,2) identity
   graph on $\{1,2,3,4\}$ accepted; corruption raises rank to 4 and is
   rejected (`controls.tensor22_*`).
2. **Duplicated-column spark** — clone keeps rank 4 at (2,4) on a 4-set
   (`controls.dup.spark4`, `controls.base.spark4`).
3. **Concat-vs-Kronecker guard** — true dimension 8 = $2\cdot4$; the
   5-dimensional block-concatenation impostor rejected
   (`controls.kron_guard.*`).
4. **Discriminating ambient plant** — a verified size-9 matching circuit
   (predicate (13): rank 8, all nine deletions rank 8) plus one extra cell
   gives rank EXACTLY 9 in the 9-dimensional ambient (`not merely ≤`)
   (`controls.nine.*`, `controls.ambient10.exact_rank`); the two-A-fiber
   nonminimal 10-set is rejected through a failing deletion
   (`controls.nonmin10.reject_deletion`).
5. **Non-Möbius plant** — an all-distinct 6-set at (2,4) GF(13) off every PGL
   graph is independent with rank 6 — square $\det\ne0$
   (`controls.nonmobius_independent`); wrong-coefficient $q^2(q-2)$ and
   selected-pole replays pass (`controls.wrong_coefficient_*`,
   `controls.selected_pole_*`).
6. **C2 witness plants** — the duplicate-tailed 9-set is dependent with rank
   < 9 (REJECT path for corank-one-with-singular-deletion); a genuine
   grid ten-set with all deletions rank 9 exists and is ACCEPTED
   (`controls.witness9.dependent`, `controls.c2_accept10_exists`).
7. **C3 plants** — verified CI 8-set circuit accepted with pencil nullity 2
   and NO common component (both Sylvester resultants nonzero);
   shared-$(1,1)$-component pencil detected and rejected
   (`controls.ci_*`, `controls.shared_component_reject`).

Runtime: `nice -n 10` asserted; `RLIMIT_CPU` (5400, 5400) set and asserted;
all five thread caps = 1; bytecode guards both active. Budget: CPU 1627.5 s
of 5400; wall 1756.1 s of 6000. Two independent rank algorithms agreed on
every computation (dual-rank cross-check inside `dual_rank`).

## 6. Defects (disclosed, remediated)

See `defect_log.json`. Summary: five pre-freeze launch aborts on instrument
bugs (a diagonal-tie rank probe; graph supports emitted as bare x-tuples; a
key-error in the p=13 want map; variable shadowing in a summary print; a
p=7 expectation pinned from a probe comparing residues to $\{1..7\}$ instead
of $F_7$ — the corrected expectations are $588$ at $n=7,k=6$ and $42$ at
$(3,4,7)$, both then re-derived independently by the PGL sum with the folded
residue domain). All broken outputs preserved. A post-run audit found the
preregistered gate-PGL replays missing from the first complete run's battery;
the complete instrument was rerun in full (second complete run, 103 asserts,
exit 0). No census was ever recorded from an abortive launch; no assertion
was weakened; the prereg was not amended.

## 7. Exact boundary of certification

Certified: Theorem T transport with obligations (a)–(c); C1 forward
dependence law and its per-prime support-set PGL equalities at all 15
registered anchor cells; C1 sub-$k$ emptiness at the registered anchors; C2
at sizes $r_Ar_B$ and $r_Ar_B+1$ with the gate-9 predicate re-expressions
(including the $5\times5$ GF(7) 36378/181960 reproduction); C3's proved
no-common-component direction with the gate-10 bootstrap counts 560/416,
2192/1344, 1632/928.

Not certified (scoped in the prereg, unchanged): C1's converse; C1 emptiness
below $r_A+r_B$ beyond the registered anchors; the C3 no-fixed-component
characterization of the whole slice; arbitrary $X,Y$ closed forms beyond the
PGL sum; unregistered spark values; three or more factors; extension fields.

## 8. Runtime

Final complete run: 103 assertions, exit 0, CPU 1627.534 s / 5400 s cap,
wall 1756.066 s / 6000 s cap, max RSS recorded in
`controls_results.json`, under `nice -n 10` with `RLIMIT_CPU` hard-limited
at 5400 s and all thread caps pinned to 1.
