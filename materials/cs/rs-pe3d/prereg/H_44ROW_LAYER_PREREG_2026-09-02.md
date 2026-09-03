# PREREG-RSPE3D-16-44ROW-LAYER — gate **H-44ROW-LAYER**

Written 2026-09-02 UTC by agent `RsPe3dH2`, immediately after
H-GRS-GENERAL-SIZE was frozen, closed FROZEN-CERTIFIED (run
`20260902T113614Z_a999d1b0_5daa598b6a04`), and independently verified by Main
at a new prime. This file is path-scoped committed before any new
parameter-dependent computation. Binding order: preregistration -> commit ->
`campaign.py init` -> byte-identical in-run copy with source commit and
SHA-256 -> controls before censuses -> analytic record -> `freeze` -> exactly
one `close --verdict`. Write only under `cs/rs-pe3d/`; do not edit any frozen
run or root ledger.

## 0. Scope and hard cap

Target: the tied four-row layer $(4\times n)\otimes(4\times n)$ GRS/Vandermonde
pairs, $d_A=d_B=5$ (factor spark $5$), ambient dimension $r_Ar_B=16$, circuit
sizes $5..17$. This gate is a **layer-by-layer prediction test of the frozen
general theorem**, not an exhaustive classification: sizes 9–17 are OUT OF
SCOPE. **Hard cap 5400 CPU seconds / 6000 wall seconds** (Main's 90 CPU-min
order), `RLIMIT_CPU` soft→hard $(5400,5400)$, `nice -n 10` asserted
(`os.nice(0)==10`), all five thread caps = 1, bytecode guards both active,
exact GF(p) integer arithmetic only.

Registered instances:

1. **(a) size 5 (fiber layer):** certify "size-5 circuits of $A\otimes B$ are
   exactly A- or B-fiber 5-sets" at $n=5$ and $n=6$ over GF(11), GF(13);
   counts $2n\binom n5$.
2. **(b) size 6 (= $d+1$) empty at $n=5,6$:** Thm-CRIT (n-factor gate) lets a
   size-$(d+1)$ tied-profile circuit exist only in the two-factor nonfiber
   channel; the crossing channel sits at $d_A+d_B-2=8\ne6$, and fixed-
   coordinate lifts require factor circuits of size 6 — none for rank-5 MDS
   factors (only size-5 circuits, gate-7/crossing taxonomy at $(5,5)$
   profiles aside, which are lifts of factor circuits and are excluded since
   a 5-row MDS Vandermonde has no size-6 circuit). Predicted
   **empty**; certify by exhaustive census at $n=5$ (all $\binom{25}{6}$
   supports) and $n=6$ (all $\binom{36}{6}$ supports) with deletion witnesses.
3. **(d) size 7 empty at $n=5,6$:** between $d+1=6$ and the crossing size 8,
   with no fiber channel (a factor 7-circuit needs spark $\le7$: for a $4\times n$ MDS factor spark is $5<7$, impossible); predicted **empty**;
   certify by exhaustive census at $n=5,6$ ($\binom{36}{7}=8.35$M — fitted to
   the cap, CPU reported).
4. **(c) size 8 (crossing + first all-distinct layer):** predicted circuits
   = crossings $25\,C_A(5)\,C_B(5)$ (i.e. $25\binom n5^2$ on an $n\times n$
   grid; note $(5,5)$-profile crossings have $|R|=|J|=5$) **plus** the
   all-distinct Möbius-8 restrictions with count
   $\sum_M\binom{k_M}{8}$ (C1 at $r_A+r_B=8$). Verification:
   - **as sets at $n=8$:** every all-distinct 8-set is a bijection $U\to V$;
     sweep all $8!=40320$ per $(U,V)$ pairing structure and compare the
     circuit set against $\sum_M\binom{k_M}{8}$ **as support sets** over two
     primes (GF(13) primary, GF(11) secondary);
   - **exhaustive tied census at $n=5$:** all $\binom{25}{8}=1{,}081{,}575$
     supports, deletions only where rank dips, to confirm NO other size-8
     profile appears beyond fibers/crossings/all-distinct (CPU reported);
   - **$n=6$ full exhaustiveness (30.3M) is scoped out under the cap.**
     Instead: construct all $25\binom 65^2=900$ crossing supports at $n=6$
     and verify each is a circuit, and verify the Möbius-8 set membership
     for all-distinct 8-sets sampled and for the full bijection sweep at
     $n=8$; the un-tested tied remainder at $n=6$ is recorded as a scoped
     non-claim in the theorem file.
5. Sizes 9–17: no claim, no sweep.

## 1. Predictions being tested (all derived from frozen gates, cited not
rederived)

- **L5 (size 5 = fibers).** From the general theorem (H-33ROW-OPEN §2 P1
  analog at $r_A=r_B=4$): ties lift only factor circuits, and a 4-row MDS
  factor has only size-5 circuits, so size-5 circuits of $A\otimes B$ are
  exactly $2n\binom n5$ fiber 5-sets. **Certify at $n=5,6$** (exhaustive
  $\binom{25}{5}$, $\binom{36}{5}$) as support-set equality.
- **L6 (size 6 empty).** $d+1=6$; Thm-CRIT requires either a
  fixed-coordinate lift (factor size-6 circuit — impossible, spark 5) or the
  two-factor tied channel with $d_A+d_B-2 = 8\ne6$. **Empty at $n=5,6$.**
- **L7 (size 7 empty).** No channel: lifts need factor size-7 circuits
  (impossible), crossing needs $8$, all-distinct needs $r_A+r_B=8$.
  **Empty at $n=5,6$.**
- **L8 (size 8 = crossings + Möbius-8).** Crossings
  $25\binom n5^2$; all-distinct restriction layer
  $\sum_M\binom{k_M}{8}$ verified as sets; no third profile (tied census at
  $n=5$).

## 2. Control battery (before any census; both directions)

1. **Pure-tensor ACCEPT/REJECT:** an A-fiber 5-set (GF(13), verified) is
   accepted by (C); a corrupted non-tensor column must raise rank and be
   rejected.
2. **Discriminating ambient plant (Main's standing order):** a 17-set
   (16-column tensor basis + 1) is dependent with rank **exactly 16** (not
   merely ≤) with a per-subset witness; a 16-set with one singular deletion
   (two swapped-fiber halves) must FAIL (C2)(i)-style minimality and be
   rejected; the 16 tensor-basis columns themselves have rank exactly 16.
3. **Non-Möbius 8-set plant:** an all-distinct 8-set off every PGL graph
   (checked against the $\binom{k_M}{8}$ enumeration) must be independent
   (rank 8) at GF(13).
4. **Concat-vs-Kronecker guard:** true dimension 16 = $4\cdot4$; a
   7-dimensional block-concat impostor rejected (coordinate check).
5. **Crossing plant:** a verified $(5,5)$ crossing support is a circuit;
   deleting its centre cell drops dependence (fail).
6. **Duplicate-column spark:** clone keeps spark; base spark 5 asserted.
7. **Möbius-8 negative control:** a proper all-distinct 7-subset of a
   verified Möbius-8 circuit must be independent (rank 7) — confirms
   minimality of the restriction layer.

## 3. Stage plan for `run_44row_layer.py`

- **A — transport:** reuse the verified Kronecker/evaluation machinery;
  dual-rank cross-checks everywhere.
- **F — control battery** (§2) **before censuses.**
- **B — PGL layer at GF(13), GF(11):** canonical enumeration
  ($p(p^2-1)$ asserted), per-map $k_M$, sums $\binom{k_M}{8}$ where $k_M\ge8$.
- **C — L5 census:** $n=5,6$ exhaustive, support-set equality with fibers.
- **D — L6/L7 census:** $n=5,6$ exhaustive emptiness (deletion witnesses on
  any rank dip; CPU reported).
- **E — L8:** $n=5$ exhaustive tied census (1.08M, CPU reported) +
  crossings-vs-rest partition; $n=8$ all-distinct bijection sweep (40320)
  with support-set equality against the Möbius-8 predictions; $n=6$ crossing
  construction (900 supports, each verified a circuit).
- **G — record,** per-stage checkpoints, budget checks per family; any
  defect freezes the family, preserves output, reruns that family.

## 4. Promotion rule and verdict

**FROZEN-CERTIFIED** iff: L5, L6, L7 hold at $n=5,6$ over both primes with
exact support-set/emptiness evidence; L8 holds at $n=5$ (no third profile)
and the $n=8$ all-distinct layer matches the PGL/graph sum as sets over two
primes; the constructed $n=6$ crossings are all circuits; and the full §2
battery passes with every planted REJECT firing. Otherwise
**FROZEN-INCONCLUSIVE**, naming the failed layer in a layer table while
certifying every layer that did pass inside the record. The $n=6$ size-8
tied exhaustiveness and all of sizes 9–17 are recorded as out of scope, never
as certified claims.

## 5. Defects, budgets, arithmetic

Exact integers mod $p$ only; inverses via `pow(x,-1,p)`; no floats/NumPy.
`sys.dont_write_bytecode=True` before all imports; launch
`nice -n 10 env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
NUMEXPR_NUM_THREADS=1 python3 run_44row_layer.py` inside the run dir.
`resource.setrlimit(resource.RLIMIT_CPU,(5400,5400))` before long stages.
Any defect is disclosed in `defect_log.json`, broken output preserved,
affected families rerun completely. No assertion weakened; no count repinned
without recomputation; the prereg is not amended.
