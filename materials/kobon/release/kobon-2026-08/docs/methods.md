# Campaign methods

This note records the implementations actually used in the 2026-08 Kobon
campaign.  It separates mathematical implications, exact witness checking,
SAT discovery, and proof-producing SAT.  Paths are relative to
`~/jinleic-workspace`.

## 1. Decision problem and soundness direction

For fixed `n` and target `T`, `math/kobon/engine.py::build_model` constructs a
CNF relaxation of:

> there are `n` distinct affine lines and at least `T` triangles supported by
> triples of those lines whose open interiors are pairwise disjoint.

The convention is broad: parallel classes, finite multipoints, selected
triangles crossed by unselected lines, and shared boundaries are allowed.
Only overlap of **open** interiors is forbidden.

Every structural clause is necessary for a real line arrangement.  Therefore a
proof-checked UNSAT result rules out a real arrangement.  The converse is not
assumed: a SAT assignment is only a combinatorial candidate and becomes a lower
bound only after exact realization and `verify_selection`.

## 2. Main SAT variables

Line labels are sorted by slope; equal slopes form contiguous parallel classes.
The principal variables are:

| family | meaning when true |
|---|---|
| `P(i,j)` | lines `i,j` cross at a finite point; false means parallel |
| `C(i,j,k)` | the three finite pairwise crossings coincide (a concurrent triple) |
| `X(r,i,j)` | on line `r`, its crossing with `i` strictly precedes its crossing with `j` |
| `S(i,j,k)` | the nondegenerate triangle supported by those three lines is selected |
| `TC(i,j,k,r)` | selected triangle `(i,j,k)` is crossed through its open interior by outside line `r` |

Do not confuse the Boolean family `C(i,j,k)` with the scalar crossing count
usually denoted `C_x` in the capacity and face inequalities.  `TC` literals sum
to that scalar line-through-triangle incidence count.

Auxiliaries make the geometry propagate in both directions: `U(i,j)` orders two
parallel lines; `SA/SB(p,q,r)` records the strict side of vertex `P_pq` relative
to line `r`; `B(r,u,a,b)` says the crossing with `u` lies strictly between the
crossings with `a,b` on `r`; and `IN(p,q,T)` flags a vertex strictly inside a
triangle.  `ZP/ZU/ZV` are Tseitin representatives used to charge each
multipoint exactly once.

### Structural clauses

The `A1`--`A7` groups enforce contiguous/transitive parallel classes,
concurrency implying crossings, concurrency closure on four labels, strict
crossing orders with concurrency ties, order transitivity and tie substitution,
the three-line order identity, and consistent ordering of parallel lines.

Selection first implies three finite pairwise crossings and no concurrent
support triple.  The `R1`--`R4` groups then impose necessary pairwise
open-interior disjointness conditions: disjoint side intervals on a shared line,
no strict vertex containment, no proper side-by-side transversal crossing, and
the concurrency boundary-contact case.  Exact `SA/SB/B` reification prevents
arbitrary auxiliary truth values from weakening these clauses.

## 3. Cardinalities, exact selection, and totalizers

`build_model(n,T)` initially encodes

```text
sum_t S(t) >= T
```

with PySAT's sequential-counter encoding.  The production cube builders then
call `add_exact_selection`, which adds

```text
sum_t S(t) <= T
```

with `CardEnc.atmost(..., encoding=EncType.totalizer)`.  This is without loss of
generality: every larger disjoint family contains a `T`-element subfamily, and
deselecting triangles only disables constraints.  Thus production instances
select exactly `T` triangles.

The bounded-face and linewise weighted lists are encoded with sequential
counters.  The fixed-arrangement max-packing prober separately uses a totalizer
for `sum candidate_triangle >= K` at each tested `K`.

## 4. Face and per-line capacity constraints

`add_triangle_crossing_indicators` exactly reifies `TC(T,r)`: it is true iff
`S(T)` is true and the three vertices of `T` include at least one point strictly
above and one strictly below `r`.

For an essential arrangement, with `Q` parallel pairs and a finite point `p`
incident with `k_p` lines, the exact bounded-face count is

```text
F_b = binom(n-1,2) - Q - sum_p binom(k_p-1,2).
```

Each selected triangle consumes a bounded face, and every outside-line/open-
interior incidence consumes another face inside that triangle.  Pairwise
interior-disjoint triangles consume disjoint face sets.  The global encoded
inequality is therefore

```text
T + sum_(T,r) TC(T,r) + Q + sum_p binom(k_p-1,2)
    <= binom(n-1,2).
```

A `k`-fold point is represented once by its three least-labelled lines and
`u=k-3` later incident lines.  The encoded weight
`1 + 2u + binom(u,2)` is exactly `binom(k-1,2)`.

With `per_line=True`, the builder also encodes for every line `r`

```text
sigma_r + gamma_r
    <= n - 2 - q_r - sum_(p on r) k_p(k_p-4)/k_p
     = n - 2 - q_r - sum_(p on r) (k_p-4),
```

where `sigma_r` counts selected sides supported by `r`, `gamma_r` counts exact
`TC(T,r)` incidences, and `q_r` counts lines parallel to `r`.  Summing over
lines yields the audited capacity theorem

```text
3|S| + C_x + 2Q + sum_p k_p(k_p-4) <= n(n-2).
```

The face inequality and capacity inequality are independent necessary budgets.
At `n=14`, combining them gives `6T + 4C_x + 5Q <= 402`, hence the unconditional
campaign bound `K_gen(14) <= 67`.  At capacity equality
`3T=n(n-2)`, the Reduction Lemma forces `C_x=Q=0`, every multipoint to have
multiplicity four, and all selected triangles to be uncrossed arrangement
faces.

## 5. Exact-rational geometry core

All certificate geometry uses `fractions.Fraction`; no floating-point predicate
is trusted.

1. `crossings(n,ms,bs)` computes every finite intersection of
   `y=m_i x+b_i` exactly and stores `None` for a parallel pair.
2. `tri_ok(Xp,t)` requires all three pairwise intersections to exist and the
   resulting support triangle not to collapse at a multipoint.
3. `tri_verts` returns the three exact vertices.
4. `hom((x,y))` clears the two rational denominators and returns an integer
   homogeneous triple `(X,Y,D)`.
5. `orient_h` evaluates the corresponding exact determinant without rational
   division.
6. `interiors_disjoint_h` applies the separating-axis theorem in both
   directions.  Weak separation is accepted, so boundary touching and shared
   edges remain legal while an open-interior overlap is rejected.

`verify_selection` validates line counts and uniqueness, triangle index shape
and uniqueness, every `tri_ok`, and every selected pair with the homogeneous
separating-axis oracle.  The Maiorana verifier adds a separate exact census and
a cevian-aware test for an outside line that passes through a triangle vertex
and then enters its open interior.

## 6. Monoliths and cubes

A **monolith** leaves every `P(i,j)` and `C(i,j,k)` free and asks the solver to
search all degeneracies admitted by the relaxation.  In
`scratch/kobon/n14/concurrency_cases.py`, the recipe name is `n14monolith`; a
regression assertion checks that the builder emits no pattern unit clauses.

A **cube** fixes an exact parallel-class pattern with `P` units and may plant
one or more positive concurrency triples.  A planted triple represents its
orbit under the symmetries that preserve the normalized parallel pattern.
Richer multipoints nest inside these cubes because any `k>=4` point makes
several triple-concurrency literals true.  Current concurrency cubes do not add
negative `C` units merely to exclude richer signatures.  Only one globally
justified crossing-order symmetry bit may be fixed; extra per-cube `X` units
are not assumed WLOG.

The strategy is:

1. use the capacity, face, and pair budgets to enumerate admissible degeneracy
   signatures;
2. normalize parallel patterns and concurrency orbits;
3. run cheap cubes for discovery and branch triage;
4. retain a true monolith as the assumption-free check;
5. promote an UNSAT branch only after its CNF and proof satisfy the discipline
   below.

A handful of built or running cubes is not automatically an exhaustive cover.
Coverage must be proved separately by the pattern/orbit ledger.  Current live
n=14, n=18, and n=20 Kissat/CaDiCaL runs without checked proofs remain
`DISCOVERY-ONLY` or `PENDING`, regardless of solver progress.

## 7. DRAT discipline

The campaign uses the following non-negotiable promotion rule:

1. Pin the exact DIMACS CNF (header, byte count, and preferably SHA-256) and the
   exact generator revision/command.
2. Run a proof-producing SAT solver.  Exit code `20` or the text
   `s UNSATISFIABLE` is **not** by itself evidence.
3. Require a complete DRAT stream; an in-flight, truncated, or killed stream is
   not a certificate even if it is tens of gigabytes long.
4. Check the exact CNF/DRAT pair with an independent `drat-trim` invocation and
   require `s VERIFIED`.
5. For a split proof, byte-check that every subcube is the parent plus the stated
   unit clauses and prove that the polarities form a tautological cover; then
   verify every subproof independently.
6. Record verifier output and exit status.  Only then may `UNSAT` be upgraded
   from discovery to machine-checked proof evidence.

This is the discipline that repaired the historical `n=10,T=26` proof family:
nine complete cube proofs were rechecked, while two truncated cubes were
re-proved by exhaustive three-literal splits and all sixteen replacement DRATs
were independently verified.

## 8. Fixed-arrangement max-packing prober

`scratch/kobon/max_packing.py` answers a different, finite question: the exact
maximum for one supplied rational arrangement.

- Vertices of the overlap graph are all triples passing `tri_ok`.
- An edge joins two triples exactly when `interiors_disjoint_h` says their open
  interiors overlap.
- A set of selected triangles is therefore an independent set.
- The SAT base has one Boolean per candidate and one binary conflict clause per
  overlap edge.
- At each `K`, a totalizer requires at least `K` vertices.  A SAT model gives a
  `K`-packing; UNSAT at `K+1` certifies optimality for that fixed arrangement.

The written `*.maxpacking.json` is conclusive only when
`"optimality_proved": true`.  The program deliberately exits normally on a
time-budget stop, so shell exit `0` alone must never be reported as an optimum.
The completed census proves equality with the stored count for all 15 Maiorana
54-witnesses and for the Bader/other `n=14` 53-witnesses and the known
`n=11,18,20` witnesses (32, 93, 116 respectively).

## 9. Direct-gap endpoint frontier

`scratch/kobon/gap_faces.py` has opt-in, monotone strengthening flags. Legacy
sector and shared-ray outputs are preserved byte-for-byte.

- `--endpoint-closure` extends the shared-ray side argument to faces sharing
  only one line at a multipoint endpoint and enforces opposite apex sides when
  two selected faces share an arrangement edge.
- `--simple-bound u` reifies `MI(t)`, meaning selected face `t` has at least
  one finite multipoint vertex, and requires at least `T-u` such faces. The
  caller must supply a proved same-\(n\) simple-arrangement upper bound.
- `--k4-bound` emits the four prime clauses for `K(4)=2` on every four-line
  subset.
- `--reify-faces` makes every `S(t)` equal the exact direct-gap predicate and
  counts all true faces. It is sound but not used in the production lane
  because \(n=9\) controls regress sharply.
- `--chirotope-gp` adds zero-aware rank-three projective determinant
  constraints. Exact audits pass, but no semantic gap was found through
  \(n=7\), and the \(n=8\) A/B run regressed; it remains propagation-only.

For \(n=12,T=39\), `--simple-bound 37` is justified by the exact simple
maximum and forces two multipoint-incident selected faces. The two current
endpoint CNFs and their promotion gates are pinned in
`scratch/kobon/frontier_endpoint_research.json`.
