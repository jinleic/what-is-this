# Tied four-row layer prediction test, sizes 5–8
## gate H-44ROW-LAYER, run `20260902T135324Z_516630ab_ed03d008f0f6`

Agent `RsPe3dH2`, 2026-09-02 UTC. Preregistration:
`prereg/H_44ROW_LAYER_PREREG_2026-09-02.md`, path-scoped source commit
`92daf3507e981b278cd48797bc91aedcd9ff6855`, SHA-256
`979b910ab7f770d5cc9e5a50b568cb95218f3b891a5f96bd4222501993daa957`. The
in-run `pre_statement.md` is byte-identical to that source (`cmp` and
SHA-256 checked); the preregistration was not amended and no
parameter-dependent compute preceded the commit.

Setting: GRS/Vandermonde pairs $(4\times n)\otimes(4\times n)$,
$d_A=d_B=5$ (factor spark 5), ambient $r_Ar_B=16$, circuit sizes 5..17.
This gate is the layer-by-layer prediction test ordered by Main: sizes
5–8 in scope, 9–17 out of scope (non-claims).

## 1. Layer table (final)

| layer | prediction | test | result | status |
|---|---|---|---|---|
| 5 | fibers only, $2n\binom n5$ | exhaustive $n=5,6$, GF(11)+GF(13) | 10 and 72 circuits = fiber set (support-set equality both primes) | **CERTIFIED** |
| 6 | empty ($d+1$; Thm-CRIT tied channel needs $d_A+d_B-2=8\ne6$; lifts need factor 6-circuits — none, spark 5) | exhaustive $n=5,6$, both primes | 0 circuits; every dependent set was 5-fiber + extra, rejected by deletion witness (200 / 2160) | **CERTIFIED** |
| 7 | empty (no channel: lifts need factor 7-circuits — impossible at spark 5; crossing sits at 8; all-distinct needs $r_A+r_B=8$) | exhaustive $n=5,6$, both primes ($\binom{36}{7}=8.35$M swept; CPU reported below) | 0 circuits; every dependent set fiber-containing and nonminimal (1900 / 31320) | **CERTIFIED** |
| 8 | crossings $25\binom n5^2$ + all-distinct Möbius-8 restrictions $\sum_M\binom{k_M}{8}$; no third profile | $n=5$ exhaustive tied census (1.08M); $n=6$ 900 constructed crossings each verified; $n=8$ bijection sweep (40320) as SETS vs PGL sum over both primes | $n=5$: 25 circuits = 25 crossings, all other dependencies fiber-nonminimal (11400); $n=6$: 900/900; $n=8$: 8 circuits = PGL sum 8 = graph set (GF(11)), 4 = 4 (GF(13)), restriction-minimality control passed | **CERTIFIED within registered scope** |
| 8 at $n\ge6$, tied exhaustiveness | out of cap (30.3M) | — | — | **SCOPED OUT (non-claim)** |
| 9–17 | — | — | — | **OUT OF SCOPE (non-claim)** |

## 2. Exact numbers

- **L5:** $n=5$: 10 = $2\cdot5\binom55$; $n=6$: 72 = $2\cdot6\binom65$; at
  GF(11) and GF(13); support-set equality against the constructed fiber
  family (`L5.*.set_equality`).
- **L6:** dependent rank-5 six-sets: 200 ($n=5$), 2160 ($n=6$) — every one a
  5-fiber + lone cell (fiber_extra = dependent), each rejected by an
  exhibited singular deletion (deleting the lone cell leaves the dependent
  fiber); 0 circuits at both primes.
- **L7:** dependent rank-6 seven-sets: 1900 ($n=5$), 31320 ($n=6$) — all
  fiber-containing, all nonminimal; 0 circuits at both primes. The
  $\binom{36}{7}$ sweeps are the cap-dominant item.
- **L8:** $n=5$ census: 1,081,575 supports; dependent 11425 of which 11400
  fiber-extra (nonminimal) and 25 = exactly the $(5,5)$ crossings
  (`L8.n5.*`, `other == circuits == 25` at both primes). $n=6$: all
  $25\binom65^2 = 900$ crossing supports constructed and each verified a
  circuit (`L8.n6.*.crossings_all_circuits`). $n=8$: all-distinct
  dependent = circuits = PGL sum = graph support set: 8 over GF(11), 4 over
  GF(13) (the maps carrying all eight residues into $Y$, counted
  $\binom{k_M}{8}=1$ each); every 7-subset restriction of a verified
  Möbius-8 circuit is independent (minimality control).

## 3. Proof obligations registered and how they were met

- **Fiber exactness (L5):** factor MDS (any 5 distinct factor columns
  dependent, any 4 independent — exercised throughout) + fixed-coordinate
  lift argument; verified as support-set equality.
- **L6/L7 emptiness:** Thm-CRIT channel condition (tied channel needs
  $d_A+d_B-2=8$), lift condition (factor circuits only of size 5 = spark),
  all-distinct channel needs $r_A+r_B=8$; certified exhaustively at the
  registered $n$'s. The exact fiber-channel condition stated in the prereg:
  a $k$-cell single-column lift is a circuit iff $k$ is a factor circuit
  size, i.e. $k=5$ — no size-6 or size-7 fiber circuits.
- **Crossings (L8):** frozen Theorem X machinery (construction + converse
  for $d_A,d_B\ge3$); construction verified circuit-by-circuit at $n=5,6$;
  exhaustiveness at $n=5$ ties certified.
- **All-distinct Möbius-8 (L8):** C1 forward at $r_A+r_B=8$ from the
  predecessor gate H-GRS-GENERAL-SIZE; certified as support-set equality
  against the canonical PGL enumeration with $|\mathrm{PGL}(2,p)|$ asserted
  (1320 / 2184).

## 4. Controls (both directions, all fired)

16-set tensor basis (4×4 subgrid) rank **exactly 16**; 17-set dependent with
rank exactly 16 (`controls.ambient17.exact_rank_16`); 15+duplicate 16-set
passes corank 1 but carries a singular basis-column deletion — rejected
(`controls.wrong16.*`); A-fiber 5-circuit accepted and its corrupted variant
rejected (`controls.fiber5.*`); crossing 8-circuit verified
(`controls.crossing.circuit`); concat-vs-Kronecker guard (16 true, 7-dim
impostor rejected); base spark 4 asserted. Non-Möbius and restriction-
minimality controls live in stage E (`L8.n8.*.restriction_independent`,
set-equality against PGL graphs).

## 5. Defects (disclosed, remediated)

`defect_log.json`: (1) a control-plant arithmetic slip (wrong16 expected
rank < 15; corrected to the intended rank-15-with-singular-deletion REJECT
witness); (2) the first full launch was killed by SIGXCPU (exit 152) at the
5400 s RLIMIT inside the GF(11) $n=8$ sweep — the all-distinct enumerator
walked $\binom{64}{8}\approx4.4\times10^9$ combination tuples instead of the
$8!=40320$ bijections; fixed by enumerating the bijections directly (both
the count sweep and the support-set comparison), pilot-verified at GF(13)
(4 circuits = PGL sum). The killed launch's console log is preserved as
`run_output_killed_sigxcpu_launch1.log`; no census was recorded from it
beyond the already-checkpointed L5/L6/L7/L8-partial results, and the complete
rerun reproduced everything from the beginning. Prereg unamended.

## 6. Runtime

Final complete run: 49 assertions, exit 0, CPU 2463.8 s / 5400 s cap, wall
2563.0 s / 6000 s cap, `nice -n 10` asserted, `RLIMIT_CPU` (5400, 5400)
set and asserted, all five thread caps = 1, bytecode guards active, dual
rank cross-checks on every decision. Sweep CPU shares: L7 $n=6$ sweeps
dominate (~1250 s per full pass); $n=5$ L8 census ~270 s per prime; $n=8$
permutation sweeps ~6 s per prime after the fix.

## 7. Exact boundary of certification

Certified: L5, L6, L7 at $n=5,6$ over GF(11)+GF(13); L8 at the registered
scope — $n=5$ tied exhaustiveness, $n=6$ crossing construction, $n=8$
all-distinct set equality over both primes.

Not certified / out of scope (recorded as non-claims): tied size-8
exhaustiveness at $n\ge6$; all of sizes 9–17; any claim about point sets
other than $\{1..n\}\bmod p$ at the registered $n$'s and primes; any
structural parameterization of the dependent-nonminimal families beyond the
fiber-containment classification verified here.
