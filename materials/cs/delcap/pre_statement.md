# delcap Pre-Statement — Gate C (finite-length claims)

Committed BEFORE any Gate C computation run. Owner: DelcapGateC. Date: 2026-08-30.
Replaces nothing; the Gate A/B pre-statement in the git history stays untouched (this file
is the working pre-statement that already served Gates A/B; this section extends it).
NOTE (league rule): an earlier agent overwrote this file's Gate A content; the committed
contract for gates A/B lives in `campaigns/` manifests. This gate-C section is appended
BEFORE the first gate-C run and is frozen at commit time.

## Papers, read first-hand today from ar5iv full text (owner reads, not REPORTED)

### P1 — arXiv:2504.20961 — Morozov & Duman (2025-04-29)
"Simple Finite-Length Achievability and Converse Bounds for the Deletion Channel and the
Insertion Channel" (Bilkent). Four extracted items:
1. **Optimization solved.** The layer-oriented converse bound (LO-CVB): for a layer
   `L = union of output-length classes {F_2^w : w in Lambda}` of the deletion channel
   `D_m^(delta)` (layers are output-length classes, equiprobable for every input),
   `M(W,eps) <= L(n,m,eps,Lambda) = (sum_{w in Lambda} tau_w)^n / ((sum_{w in Lambda} p_w)^n - eps)`
   with `tau_w = E(m,w) * delta^(m-w) (1-delta)^w`, `p_w = C(m,w) * delta^(m-w) (1-delta)^w
  `, and `E(m,w) = sum_y max_x binom(x|y)` the sum of maximal embedding numbers
   (eqs 25-30). Minimize over subsets Lambda of {0..m} of output lengths with
   E(m,w) available, then over 20<=m<=32 with side information trick N=mn.
   The achievability side is a *greedy* ML-code-selection algorithm (Alg. 1, GAVB),
   exponential, run only at n=1. NOT a convex program; the LO-CVB is a closed-form
   evaluation given integer tables E(m,w).
2. **Truncation.** Complete E(m,w) tables for m<=23 (their Table I) and partial
   (w in [0..w0] union [w1..m]) for 24<=m<=32 (their Table II); insertion side
   E1(w,m) complete for m<=16, partial to 25. Their table Entries are exact integers
   claimed computed by exhaustive/DP enumeration with reverse/inverse symmetries.
3. **Stated numerical precision.** A general BA-convergence tolerance a=0.005 for
   reported C_{n,k} bounds; a=0.05 for C_{31,k}, 14<=k<=18 (paper §6). LO-CVB
   itself needs no tolerance beyond the integer tables; the code rate values in
   their Table III are printed to 5 and 6 decimal places.
4. **Where the constants appear.** Table III (delta=0.2, eps=0.2, code-rate
   LO-CVB vs BEC), Fig. 2 (delta in {0.05,0.2,0.5,0.8}, LO-CVB curves vs BEC and
   normal approx), Fig. 1/3 (optimal codes and GAVB comparisons), Table IV/V
   (insertion-channel E1 tables).

### P2 — arXiv:2604.05867 — Pinto & Ribeiro (2026-04-07)
"Improved Capacity Upper Bounds for the Deletion Channel using a Parallelized
Blahut-Arimoto Algorithm" (IST Lisboa). Four extracted items:
1. **Optimization solved.** Standard Blahut-Arimoto BA upper bound on the *exact*
   deletion channel capacity C_{n,k} = sup I(X^n;Y) (their §2.3/Alg. 1), run on a
   CUDA GPU with DP-subsequence enumeration; the conversion to BDC capacity upper
   bound C(d) <= (1/n) C_n(d) = (1/n) sum_k binom(n,k) d^(n-k)(1-d)^k C_{n,k}
   (their Lemma 2, credited to Fertonani-Duman) is a finite convex combination.
2. **Truncation.** C_{n,k} computed for all k<=n<=29, plus C_{31,k} for k<=18;
   for C_{31,k}, k>18 they used the RC Lemma-10 split recursion at s=2
   (their Lemma 10) from the C_{29,*} table (their Table 2). No output-window
   truncation: exact deletion channels only (never the union-length channel).
3. **Stated numerical precision.** BA tolerance a=0.005 for almost all
   (n,k), a=0.05 for C_{31,k}, 14<=k<=18. Reported values = BA rate + tolerance
   (an additive round-UP rule, so published numbers are valid upper bounds IF their
   stopping-criterion claim holds). GPU: RTX 5070 Ti, up to 1100 iterations at n=29,
   2500 s/iteration at k=18.
4. **Where the constants appear.** Their Table 1 (C_{29,k} upper bounds k=1..29),
   Table 2 (C_{31,k}, k=1..31), Table 3 (combined C(d) upper bounds and previous
   best [[29]] comparison; the 0.3578 at d=0.64 is the RATIO C(0.64)/0.36
   propagated by the Rahmati-Duman RD15 monotonicity to all d>=0.64).

### P3 — arXiv:2607.19559 — Tavakoli, Nguyen & Bose (2026-07-21)
"Combinatorial Capacity Bounds for the q-ary Deletion Channel" (Oregon State).
Four extracted items:
1. **Optimization solved.** A BA convergent-iterates evaluation of the finite
   per-symbol capacity C_{q,n} of the q-ary deletion channel, sandwiched by:
   lower bound  C_{q,n} >= (1-d) log2 q - h2(d)                    (their Thm 1, LB1)
   refined      C_{q,n} >= LB+ = (1-d)log2 q + (1/n)H_Bin(n,1-d) - h2(d)
                       + Δ_n(d)/n                                  (their Cor. 1)
   with  Δ_n(d) = sum_{k=1}^{n-1} w_k Φ_{k,n},
         w_k    = binom(n,k) d^(n-k)(1-d)^k,
         Φ_{k,n}= (1/(q^n binom(n,k))) sum_{x,y} N_n(x,y) log2 N_n(x,y)
   and upper bound  C_{q,n} <= (1-d) log2 q  (their Thm 1). This sandwich is the
   "d=1/2 sandwich" of the ticket: the sandwich is valid for all d but their
   Table-1 NUMBERS are at d in {0.05, 0.10, 0.20}.
2. **Truncation.** No truncation: Φ and Δ are finite sums over the whole
   (q^n x 2^{≤n}) oracle-space via the pattern-count recursion; H_Bin is an
   n=…,full-range sum. The fine structure uses only exact integers N_n(x,y).
3. **Stated numerical precision.** BA values `BAC_(q,n)` computed with the
   Blahut-Arimoto algorithm; no explicit iteration/dps stated. Printed to 3 d.p.
   in their Table I (q=2, n in {3,5,10}, d in {0.05,0.10,0.20}). They report the
   exact uniform-input rate LB+ (their Cor. 1) as an un-toleranced closed form
   (no tolerance needed: it is exact given the integer patterns).
4. **Where the constants appear.** Their Table I, all 18 rows: q in {2,3},
   n in {3,5,10}, d in {0.05, 0.10, 0.20}.

## What Gate C actually re-derives + extends (the committed grid)

### (A) Tavakoli et al. (P3) — the SANDWICH, at their (q,n,d) grid + extension
   For each (q, n, d) in their own Table-I 18 points (q in {2,3}; n in {3,5,10};
   d in {0.05,0.10,0.20} as exact RATIONALS — 1/20, 1/10, 1/5),
   recompute in certified Arb:
     LB1(q,n,d)      [Trivial sandwich lower]
     LB+(q,n,d)      [Tightened sandwich lower, from Φ_{k,n} exact integer patterns]
     UB(q,n,d)       [(1-d) log2 q]
     Cbar_{q,n,d}    [true finite-block capacity = sup over q^n inputs] via the
                     certified-BA pipeline (mpmath BA locate -> exact-rational snap ->
                     outward-rounded Arb primal+dual) on the q-ary deletion channel
   Δ_n(d) is computed EXACTLY over the pattern-count integer array N_n:
   Φ_{k,n} = (1/(q^n binom(n,k))) * sum N_n log2 N_n  is an exact rational times
   log2 of integers; evaluated in Arb outward-rounded.
   Their (q,n,d) points, plus my extension to q=4 and q=3,n=7,8 and d=1/2
   (d=1/2 = the ticket's named target for the "improvement on their d=1/2 sandwich").
   For d=1/2 the Tavakoli sandwich is (1/2)log2 q - 1  <=  C_{q,n}(1/2)  <=  (1/2)log2 q
   (with h2(1/2)=1). The certified-BA interval gives the cert-true value inside.
   Target: a certified improvement on the d=1/2 sandwich means an Arb interval fully
   contained in [sandwich-LB, sandwich-UB] and strictly narrower.

### (B) Pinto-Ribeiro (P2) — C_{n,k} upper bounds, their points + extension
   Their Table 2 rows at n=8, k in {1..8} (small n where their BA tolerance 0.005
   is the WEAKEST — their C_{8,k} values are meaningful to compare precisely
   because they enter the C(d) combined bound); extend to (n,k) they did not
   publish at n<=10 by a full certified-BA interval [I(p*), dual].
   NOTE: their published C_{n,k} values at n=29/31 need BA on 2^29x2^k; out of
   reach on this workstation (this is documented in the gate-A runway note).
   The certifiable subset committed here: (n,k) with k<=n<=12 (small-block where
   our exact dense pipeline is exact and fast); specifically (n,k) rows
   (1,1),(2,2),(3,1..3),(4,1..4),(5,1..5),(6,1..6),(7,1..7),(8,1..8) — these
   OVERLAP with their Table 2 at the k where 29<=n<=31 is their data (n=29
   table, k=1,2  Rows 1,2 do not include small n). Given that, the honest
   comparison set is: their published SMALL-n used in Lemma-2 combinations is
   empty — Pinto-Ribeiro only publish C_{n,k} for n in {29,31}. So the Gate-C
   certified rows for P2 will be the (n,k) in {5..12} numbered rows reported as
   **NEW certified values, no published analog** (a different value-add than
   gate B's Table II), PLUS the paper-level small-n trace value implied by
   Lemma 2 with the certified intervals used in the combination. At d=0.64 the
   direct comparison IS possible: recomputing C(d) <= (1/n)·(Σ_k binom·…·C_{n,k})
   at d=0.64 using the PUBLISHED Table-2/3 constants gives their 0.1288 row →
   our certified rec uses our certified Cbar_{12,k} with the Eq.-2 recursion BUT
   at weak n it would not match their value — we state this as a partial (the
   rows we can certify exactly are limited to n<=12; their n=29/31 GPU result is
   NOT re-derivable on this workstation). This is stated up-front as a
   structural limit, not a silent omission.

### (C) Morozov-Duman (P1) — LO-CVB integer table rows + a new certified route
   Re-derive rows of their LO-CVB from the exact integer equations (24)-(30):
   the LO-CVB value is a closed-form rational expression in E(m,w) tables which
   we recompute from scratch by our own exact integer embedding-count code
   (independent of their published tables) — full E(m,w) for m<=10, then bound
   the LO-CVB in Arb. This produces a *certified interval* around their point
   value (which they print to 5 decimals) at their own targets (delta=0.2,
   eps=0.2, code rate; the F_2^w layer union Max over Lambda).
   Note: their tables only reach m<=32 (deletion); our machine reaches m<=~26 in
   the dense regime; the commit here is m in {5, 10, 15, 23} with their exact
   Table I rows (5,2)=32, (5,3)=52, (5,4)=54 — values duplicated in their
   eq (28) row that name a specific E(m,w) for m=5 — used as the ANCHOR.

## Precision & rounding policy (fixed BEFORE run)

- Arb outward-rounded balls (python-flint 0.9.0), working precision 400 bits
  per certificate (5e-140 relative-bounds); all exact rationals fmpq.
- mpmath 150-dps for the BA LOCATING step only. No mpmath number is
  used in any verdict.
- Every certified claim carries its interval WIDTH as a first-class number.
- Target interval width: absolute <= 1e-6 on C_{q,n,d} per symbol (achievable
  for n<=12 by this pipeline; the mpmath-BA dual gap is ~1e-11 n<=4, ~3e-8 n=8,
  and the snap denominator 2^40 will dominate). Wider mid-n cases expected
  ~1e-7 and reported AS MEASURED, not clamped.
- No interval iteration anywhere; the only interval arithmetic is one-shot
  outward-rounded evaluation of an explicit rational formula.

## Stop rule (pre-registered)

  - REPRODUCE / CONFIRM  : paper printed value lies inside the certified interval
                           (for P2/P1, allowing for their additive tolerance a>=0).
  - DISAGREE             : certified interval EXCLUDES the printed value.
        -> freeze evidence, `hub`-message Main, STOP that row. No README write.
  - IMPROVE              : certified interval strictly narrower than the paper's
                           (P3 d=1/2 target: interval inside the sandwich AND
                           width < sandwich width).

## What will NOT be swept (rule 7)

- Pinto-Ribeiro n in {29,31}: GPU pipeline, out of scope on this workstation.
- Morozov-Duman m>32 (deletion table), m>25 (insertion table).
- Tavakoli n>10 for the STRUCTURAL Δ_n: Δ_n needs (q^n x q^k) full enumeration;
  at q=2, n=10 it is already 1024 x 2^10 and full pattern table is fine; at
  n=12 it is 4096 x 4096 - still fine, that IS in scope; n>12 dropped.
- Any claim about the ASYMPTOTIC deletion channel capacity: closed, except
  through P2's RD15 monotonicity chain at the paper's own d=0.64 point,
  which is out of our reachable certified-n regime.

Signed: DelcapGateC. No Gate C computation run yet as of commit time.
