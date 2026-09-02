# SCOPE NOTE — Theorem X after H-MIX-H2-CLOSURE
# original run `20260901T132343Z_4f169cbb_2f81dd8f2d8d`;
# H2 closure run `20260901T135534Z_5678c46e_232421125d33`
# all-distinct successor run `20260902T011204Z_7ece99fe_a108b246f87d`
# PGL-count successor run `20260902T012629Z_2eb59adb_0aa7425efe1b`
# n-factor d+1 successor run `20260902T014115Z_25757c71_9c6065a4173c`
# two-factor d+2 successor run `20260902T021721Z_7629001e_215016bf40c3`
# two-row capstone / 2-by-3-row frontier run `20260902T023120Z_8ac826db_96f766bfe6cc`
# 2-by-3-row determinantal closure run `20260902T025427Z_1b1f202c_6247a0994bce`
# tied 3-row-by-3-row open run `20260902T031527Z_2d19c5c4_718416983678`
# tied 3-row all-distinct size-8 residual run `20260902T041030Z_64e0c4ef_b5de16e49f8a`

Status: the first six gates, the eighth gate, and the tenth gate are
FROZEN-CERTIFIED; the seventh exploratory gate and ninth tied-three-row gate
remain historically FROZEN-INCONCLUSIVE (2026-09-01/02 UTC). All ten runs are
frozen and untouched. The eighth superseded the seventh gate's 2-by-3-row
determinant gaps. The ninth closed tied-three-row sizes four through seven and
recorded exact predicates through the ambient cap; the tenth now explains its
registered all-distinct size-eight residual as the reduced complete-
intersection channel. Repeated-index size-eight templates and structural
counts above size eight remain outside that final gate. This target-level note
supersedes prior open-frontier caveats without rewriting any frozen theorem
file.

## Proved in general (all factor pairs over any odd GF(p), no swept-config
dependence)

1. **Forward construction.** For $R \in \mathrm{Circ}_A(d_A)$, $i_0 \in R$,
   $J \in \mathrm{Circ}_B(d_B)$, $j_0 \in J$, the matrix $\Gamma = c\,e_{j_0}^{\mathsf T}
   + e_{i_0}\,\delta^{\mathsf T}$ (with $\delta$ scaled so $c_{i_0} + \delta_{j_0} = 0$)
   lies in $\ker(A \otimes B)$ and has support exactly
   $S(R,i_0,J,j_0) = (R \setminus \{i_0\}) \times \{j_0\} \;\cup\; \{i_0\} \times
   (J \setminus \{j_0\})$.
2. **Circuitness + size law + relation space.** $S(R,i_0,J,j_0)$ is a genuine
   circuit of size $d_A + d_B - 2$ with profile $(d_A, d_B)$; its relation
   space is 1-dimensional, spanned by $\Gamma$ (rank proof: the two arm spaces
   $W_A \otimes b_{j_0}$ and $a_{i_0} \otimes W_B$ intersect exactly in the
   centre line; relation-space proof: dimension $= |S| - \mathrm{rank}$ and
   full support of $\Gamma$).
3. **Closed form, via injectivity of the parameterization.**
   $N_{\mathrm{mix}} = d_A\, d_B\, C_A(d_A)\, C_B(d_B)$ holds whenever the
   map $(R,i_0,J,j_0) \mapsto S(R,i_0,J,j_0)$ is injective — proved: for
   $d_A, d_B \ge 3$ (multi-cell column/row identification) and for exactly one
   of $d_A, d_B = 2$; **exactly 2-to-1** at $(2,2)$, giving
   $N_{\mathrm{mix}} = 2\,C_A(2)\,C_B(2) = \tfrac12 d_A d_B C_A(2) C_B(2)$ there.
4. **Field independence for fixed factor circuit spectra.** The mixed count is
   a function of $(C_A(d_A), C_B(d_B))$ only; if two fields give the same
   factor circuit spectra, they give the same mixed count (construction and
   injectivity are field-free given full-support circuit relations).
5. **Reach/coherence.** The channel emits size-$(d+1)$ circuits iff
   $d_A + d_B = d + 3$, saturating run-4's Thm-CRIT necessary condition.
6. **Lemmas P, Q, R.** Pure-piece decompositions are impossible for all
   $d_A, d_B \ge 2$; $(1,k)$ and $(k,1)$ shapes ($k \ge 2$) are impossible for
   $d_A, d_B \ge 3$; the $(k_1,k_2) = (1,1)$ case is exactly a crossing
   support.
7. **Converse for all $d_A,d_B\ge3$ (H2 discharged).** The successor proof
   closes the only residual $(d_A,d_B)=(3,3)$, $(k_1,k_2)=(2,2)$ corner:
   overlap is confined to the four active-row/active-column cells; tight
   support counting forces all four cells into both summand supports and
   outside the circuit support, hence they cancel; equality makes every
   component support a 3-circuit; profile $(3,3)$ forces the two third row
   indices and the two third column indices to coincide. The resulting support
   is exactly $S(T,w,Z,x)$, and its relation is proportional to the crossing
   relation. Thus P + Q + Q2 + R + H2 give the full converse for
   $d_A,d_B\ge3$, with no swept-configuration hypothesis.


## Converse status and remaining boundaries

- **Converse, main clause — now proved.** For all $d_A,d_B\ge3$, every
  size-$m_*=d_A+d_B-2$ profile-$(d_A,d_B)$ circuit is a crossing support.
  H2 is no longer an open question. The analytic closure is in successor run
  `theorem_h2_closure.md`, Steps 1–5, including the inside-$S$ overlap case,
  the equality analysis $P\cap Q=X$, the exclusions
  $w_i\notin\{u_1,u_2\}$ and $x_i\notin\{v_1,v_2\}$, and the spark-to-circuit
  implication. Exact in-run corroboration swept every normalized $(2,2)$
  parameter choice over GF(7), GF(11), and GF(13) for two distinct factor
  pairs; every profile-$(3,3)$ support was a member of the constructed
  crossing family (set equality: 144 and 90 supports per configuration, zero
  new supports), and the prior GF(13) 1,820-subset census was reproduced with
  144 measured profile-$(3,3)$ circuits equal to 144 constructed supports.
- **Spark-2 boundary.** The branches $d_A = 2$ and/or $d_B = 2$ of the
  converse (e.g. permutation-shaped sets when one side has spark 2) are OUT
  OF SCOPE of this gate (prereg §6): not proved, not measured. The 2-to-1
  counting clause at $(2,2)$ is proved as a parameterization statement only.
- **All-distinct 2-row GRS channel — classified and PGL-counted.** Gate
  H-ALLDISTINCT-CROSSRATIO proves that a size-four profile-$(4,4)$ support for
  columns $(1,x)\otimes(1,y)$ is a circuit iff
  $\det[1,y,x,xy]=0$, iff the paired points lie on one nondegenerate
  bilinear/Möbius graph, iff the pairing preserves cross-ratio. Successor
  H-ALLDISTINCT-PGLCOUNT proves the exact arbitrary-finite-set principle
  $$N_{\rm all}(X,Y)=\sum_{M\in{\rm PGL}(2,p)}
  \binom{|X\cap M^{-1}(Y)|}{4}.$$
  It derives closed forms for $X=Y=F_p$ and $X=Y=F_p^*$ by pole/zero classes;
  three-route product/CR/PGL support equality reproduces full-field counts
  200/5880 at $p=5,7$ and multiplicative counts 8/1080/117600 at
  $p=5,7,11$. Field dependence is exactly the PGL incidence profile of the
  point sets. Still open: further symbolic simplification for arbitrary
  unstructured $X,Y$ and any higher-row analogue; when $r_Ar_B>d+1$ the safe
  proven statement remains a bi-Vandermonde rank condition, not one
  determinant.
- **Three or more factors at $|S|=d+1$ — now closed.** Successor gate
  H-DP1-THREEFACTOR proves for any finite number of factors with nonzero
  columns and all $d_i\ge3$ that every size-$(d+1)$ circuit varies in at most
  two coordinates. The missing two-factor repeated-mixed profiles are
  excluded analytically by a projection floor and the repeated-edge rank sum
  $\rho_A+\rho_B\le d+1$. Hence every circuit is uniquely a fixed-coordinate
  lift of a size-$(d+1)$ factor circuit or a nonfiber two-factor circuit, and
  its exact count is
  $$\sum_i\Bigl(\prod_{k\ne i}s_k\Bigr)C_i(d+1)
  +\sum_{i<j}\Bigl(\prod_{k\notin\{i,j\}}s_k\Bigr)N^{\rm nf}_{ij}(d+1).$$
  Exact three-factor support-set equalities were 81 for $V_3^{\otimes3}$,
  1824 for $V_4^{\otimes3}$, 624 for $M_{34}\otimes V_4\otimes V_4$, and
  156 for $M_{34}\otimes V_3\otimes V_4$, always zero genuinely 3-D.
  General $n\ge3$ classifications above $d+1$ remain open; the special
  two-factor 2-row $d+2$ layer is closed next.
- **Two 2-row GRS factors at $|S|=d+2=5$ — now closed.** Gate
  H-DP2-SIZE5 classifies every five-circuit by nine bipartite path/cycle
  templates plus zero, one, two, or five all-distinct 4-minor tests. The tight
  $(1,2)/(2,1)$ channels give field-free counts
  $N_{33}=9\binom n3^2$,
  $N_{34}=N_{43}=72\binom n3\binom n4$, and
  $N_{35}=N_{53}=18n\binom n3\binom{n-1}4$.
  The remaining exact predicates include
  $N_{44}=12(M_4-N_{\rm all}^{(4)})$ and two-/five-minor sums for
  profiles $(4,5)/(5,4)/(5,5)$. The correct general decomposition support law
  is $|S|=T_1+T_2-r-c$; the simpler `-2c` expression applies only when all
  overlaps cancel and all active relations have minimum support. Full
  direct/structural/constructed support-set equality held for $n=4,5$ at
  $p=7,11,13$ (size-five totals for $n=5$: 17124/17560/17880).
  Size six and above remains open.
- **Full two-2-row capstone proved; 2-by-3-row frontier opened but not
  closed.** Successor gate H-2ROW-COMPLETE-3ROW-OPEN combines the frozen
  layers with the ambient-rank lemma $|S|\le r_A r_B+1$. Thus two 2-row GRS
  factors have exactly size-3 fibers, size-4 crossings/all-distinct
  cross-ratio supports, and the complete size-5 graph/minor families; no
  circuit has size at least 6, and size 5 attains the general bound. It also
  records the correct unequal degree convention
  $\deg_A\le d_B-1,\deg_B\le d_A-1$ outside the minimum-fiber exception.
  For $2\times n$ by $3\times n$ Vandermonde factors, exact GF(7)/GF(13)
  censuses for $n=4,5$ prove the size-5 crossing count
  $12\binom n3\binom n4$, three field-free size-6 constructions with counts
  $18\binom n3\binom n4$, $144\binom n4^2$, and
  $180\binom n5\binom n4$, and the size-7 $K_{2,3}+K_2$ count
  $12\binom n3\binom n4$. The whole gate is FROZEN-INCONCLUSIVE because the
  size-5 profile-$(5,5)$ determinant channel, size-6 profiles
  $(4,5)/(5,5)$, and residual size-7 channels lack general analytic
  parameterizations. The two-row capstone remains a proved certified
  subresult of that frozen run.
- **The 2-row-by-3-row determinantal residual is now closed.** Successor gate
  H-3ROW-DETERMINANTAL proves that a size-5 profile-$(5,5)$ matching is a
  circuit exactly when it is the graph of a unique Möbius map, with count
  $\sum_M\binom{k_M}{5}$. For six cells with one repeated B value,
  $D_6$ factors into nonzero difference/Vandermonde terms times
  $D_4=\det[1,y,x,xy]$ on the four singleton edges. This yields
  $N_{45}^{(6)}=6(|Y|-4)\sum_M\binom{k_M}{4}$ and
  $N_{55}^{(6)}=4\sum_M\binom{k_M}{4}
  ((|X|-4)(|Y|-4)-(k_M-4))$, including the necessary central-deletion
  exclusion in the latter. Every size-7 channel now has a
  necessary-and-sufficient list of nonzero six-deletion minors; the $q=4$
  channel is field-free with
  $N_{p4}^{(7)}=4p\binom{|X|}{p}\binom{|Y|}{4}g_{p-1}$.
  Direct-versus-predicted support sets agreed overall, by profile, and for all
  17 measured size-7 signatures at $n=4,5$ and $p=7,11,13$. Further
  simplification of the exact $q\ge5$ determinant sums remains open, but the
  circuit predicate is complete.
- **Tied 3-row-by-3-row factors — sizes 4 through 7 are now closed.** Gate
  H-33ROW-OPEN proves for columns
  $(1,x,x^2)\otimes(1,y,y^2)$ that size four consists exactly of fibers, with
  $N_4=2n\binom n4$ on an $n\times n$ grid; size five is empty; and every
  size-six circuit is either a crossing or the graph of a unique Möbius map.
  Thus
  $$N_6=16\binom{|X|}{4}\binom{|Y|}{4}
       +\sum_{M\in\mathrm{PGL}(2,F)}\binom{k_M}{6}.$$
  A Vandermonde-kernel/conic argument excludes every repeated and all-distinct
  size-seven profile, so the entire size-seven layer is empty. The same gate
  proves that eight points on a reduced rational degree-two graph form a
  circuit, uniquely parameterized from five pairs, with incidence count
  $\sum_{\deg R=2}\binom{|D_R(X,Y)|}{8}$; it does not exhaust size eight.
  Size nine is exactly common bidegree-$(2,2)$ curve plus all eight-deletion
  ranks equal to eight, and size ten is exactly nonvanishing of all ten
  nine-deletion determinants. Exhaustive grids covered $n=4$, sizes 4--10,
  and $n=5$, sizes 4--8, over GF(7), GF(11), GF(13), plus $n=5$ sizes 9--10
  over GF(7). Dedicated $8!$ matching sweeps left 560/416 size-eight circuits
  outside the rational-degree-two family at GF(11)/GF(13), so the gate is
  FROZEN-INCONCLUSIVE and the first structural frontier is that size-eight
  residual.
- **Tied 3-row all-distinct size-eight residual — now closed as complete
  intersections.** Final gate H-33ROW-SIZE8-RESIDUAL proves that every
  all-distinct eight-circuit lies either on a reduced irreducible
  degree-two graph of class $(1,2)$ or $(2,1)$, or is the full reduced base
  locus of a unique fixed-component-free pencil of $(2,2)$ curves. Conversely
  every support in those classes is a circuit. The graph incidence sums must
  require injectivity on the selected eight parameters; the complete-
  intersection term is the exact number of such pencils, with no unproved
  automorphism quotient. For the registered $X=Y=\{1,\ldots,8\}$ sets, both
  graph orientations are empty: every prior residual is a complete
  intersection, exactly 560 over GF(11) and 416 over GF(13). Bare rank seven
  overaccepts: the full pencil populations are 2192/1344, and the remaining
  1632/928 supports have a common $(1,1)$ component and a six-point Möbius
  subcircuit. On all 976 circuits, both elimination resultants equal the
  normalized degree-eight projection-node polynomials. There is no ninth base
  point because $(2,2)\cdot(2,2)=8$. The gate is FROZEN-CERTIFIED for this
  all-distinct channel; it does not silently classify repeated-index
  size-eight templates. The target rests for the remainder of this session.
- **Swept-config certification.** At every configuration actually enumerated
  in-run (see run `controls_results.json`, T2/T3b/T4/T5/T6/T9 set-equality
  blocks and T7 prime sweep), the machinery verified: measured circuits of
  the stated size/profile EQUAL the constructed family exactly — so all
  conclusions of §2 of the theorem file hold verbatim at those configs.

## Reader guidance

Read the predecessor theorem files together with the two records in the
seventh run, the determinantal theorem in the eighth, and the two tied-three-
row records in the ninth and tenth: original `theorem_crossing.md`, H2
successor `theorem_h2_closure.md`, cross-ratio successor
`theorem_alldistinct_crossratio.md`, PGL-count successor
`theorem_alldistinct_pglcount.md`, n-factor successor
`theorem_dp1_threefactor.md`, size-five successor
`theorem_dp2_size5.md`, seventh-run `theorem_2row_complete.md` and
`frontier_2x3row_open.md`, eighth-run `theorem_3row_determinantal.md`,
ninth-run `theorem_33row_open.md`, and tenth-run
`theorem_size8_residual.md`. The original H2/open-frontier wording is
historical: for $d_A,d_B\ge3$, Theorem X has an unconditional crossing
converse; the two-2-row classification is complete through its ambient cap;
the 2-row-by-3-row GRS residual through size seven has exact graph/minor/PGL
predicates; tied three-row products are structurally complete through size
seven with exact decision predicates through size ten; and their all-distinct
size-eight channel is exactly degree-two graphs plus reduced $(2,2)$ complete
intersections. Still open are the spark-2 converse branches, closed
simplifications of the unequal 2-by-3-row field-dependent size-seven
determinant sums, repeated-index tied-three-row size-eight templates,
structural parameterizations/counts for tied-three-row sizes nine and ten,
and general higher-rank or higher-factor classifications.
