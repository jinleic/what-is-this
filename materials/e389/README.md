# Erdős Problem #389: consecutive-product divisibility

**Canonical target.** Prove or disprove that, for every integer $n\ge1$, some
$k\ge1$ satisfies

$$
 n(n+1)\cdots(n+k-1)\mid(n+k)(n+k+1)\cdots(n+2k-1).
$$

Set $m=n-1$. The exact working form is

$$
\boxed{\binom{m+k}{m}\mid\binom{m+2k}{k}}. \tag{1}
$$

The universal statement remains **OPEN**. This directory contains exact
reductions and bounded computation, not a claimed solution.

## Current status

Labels are literal.

- **PROVED (elementary, below):** the valuation/carry and prime-power-zone
  criteria; exact band/local small-prime controls; compensation-good integer
  separation; the short-cofactor prime-power spike criterion; the guaranteed
  first-prefix repair cone and its exact size; retained-prime product
  displacement; zero-carry CRT eviction of prescribed blockers; the
  shift-spike/first-term identity; **half-window translation eviction with its
  sharp minimal offset**; and **the powersmooth necessity bound for every
  compensated prime power**.
- **PROVED (using the prime number theorem):** the short-cofactor density is
  $\log 2$, so for every fixed $m$ the witnesses have upper density at most
  $1-\log 2=0.30685\ldots$. This is a ceiling, not an obstruction.
- **PROVED (using Dirichlet's theorem):** every fixed CRT class contains
  arbitrarily large single-level bad-window obstructions.
- **CERTIFIED FINITE:** all 27 values currently listed in OEIS A375071
  ($0\le m\le26$) are witnesses. Legendre valuations and Kummer carries agree.
  Every nontrivial listed witness has minimum valuation slack zero and is tight
  at 2--39 primes. The verifier independently exhausts minimality only for
  $m=3$ ($k=207$) and $m=5$ ($k=2475$).
- **CERTIFIED FINITE:** `atlas.py` exhausts
  $1\le m\le50$, $1\le k\le50{,}000$. It reproduces the published minima for
  $m=1,\ldots,10$; no witness occurs in this rectangle for $m=11,\ldots,50$.
  This says nothing about $k>50{,}000$.
- **CERTIFIED REPAIR:** atlas schema 2 includes prime $p$ in every carry-pattern
  identity. The preserved schema-1 artifact merged bases incorrectly; an
  independent re-scan identifies every affected event and proves that witness
  truth values, per-$m$ counts, and offset-screen results were unchanged.
- **CERTIFIED FINITE:** the first 200,000,000 forward offsets from the
  published $m=26$ natural-shift target contain no $m=27$ witness. Twelve
  offsets survive the single-level sieve; complete exact checks reject all.
- **SUPERSEDED:** the adaptive eviction progression reaches a uniform
  small-prime-safe class at \(10\le t\le706\) with 32--81-digit candidates, all
  beyond the retained factorization bound. Those sizes came from demanding the
  *sufficient* uniform classes, not from eviction. The half-window translation
  evicts every old blocker at 13 digits, satisfies the exact small-prime tier
  for all 12 survivor systems with no search, and classifies all 12 completely:
  each fails only at 9--15 newly entering large primes.
- **CERTIFIED FINITE:** combining the atlas with an exact compensation-run
  sieve proves that an $m=27$ witness, if one exists, has
  $k>1{,}000{,}050{,}000$.
- **CERTIFIED FINITE:** four new contiguous exact compensation-run shards extend
  the neighbourhood of the published $m=26$ natural-shift target to the full
  interval $5{,}048{,}891{,}644{,}619\le k<5{,}049{,}891{,}644{,}619$. All
  $10^9$ candidates are rejected, $8\cdot10^8$ of them by the exact
  compensation criterion with no offset reaching the small-prime tier, and the
  longest compensation-good runs are 8, 9, 9, 8 against the required 14.
- **MEASURED EXACTLY:** the $27$-compensation-good density is $0.133708$ over
  $[50{,}015,\,100{,}050{,}027]$, $0.129223$ over
  $[100{,}050{,}015,\,1{,}000{,}050{,}027]$, and $0.123664$ in a fresh block at
  $5\cdot10^{12}$; the exact small-prime tier passes for 19,920 of 20,000
  consecutive $k$ there. Observed run spectra match a term-independence model
  within one unit of the longest run.
- **NOT PROVED (extrapolation):** regressing $\log k_m$ on $\lceil m/2\rceil$
  over the 26 published witnesses gives $R^2=0.989$ and predicts a least
  $m=27$ witness near $2.8\cdot10^{13}$; the independence model predicts
  $5.1\cdot10^{12}$. Both exceed every exhausted range here.
- **OPEN:** a uniform witness construction, a lower-bound density theorem, an
  effective global bound, or a finite shift-repair theorem.

## Exact reductions

Define the $p$-adic slack

$$
 s_p(m,k)=v_p\binom{m+2k}{k}-v_p\binom{m+k}{m}.
$$

Only primes dividing the left binomial need checking. Thus (1) is equivalent to

$$
 s_p(m,k)\ge0\qquad\text{for every prime }p. \tag{2}
$$

Kummer's theorem identifies the two valuations with base-$p$ carry counts:

$$
 v_p\binom{m+k}{m}=\operatorname{carries}_p(k,m),\qquad
 v_p\binom{m+2k}{k}=\operatorname{carries}_p(k,m+k).
$$

Consequently (2) is exactly

$$
\operatorname{carries}_p(k,m+k)\ge
\operatorname{carries}_p(k,m)\quad\text{for every prime }p. \tag{3}
$$

The checker factors only the $m$ integers $k+1,\ldots,k+m$ and $m!$ to obtain
the prime support of the left binomial. It never constructs a product or
central binomial having $k$ terms.

The retained trial-factor table is explicitly capped at inputs $10^{14}$,
which bounds its sieve at $10^7$ and covers the present certificates. This is
an implementation/resource limit, not a mathematical restriction; larger
candidates require a different exact factorizer before they are accepted.


## Natural shift framework

Let

$$
 Q(m,k)=\frac{\binom{m+2k}{k}}{\binom{m+k}{m}}.
$$

For $k\ge2$, the shift $(m,k)\mapsto(m+1,k-1)$ preserves $m+k$ and satisfies

$$
\frac{Q(m+1,k-1)}{Q(m,k)}
 =\frac{m+1}{m+2k}. \tag{4}
$$

Therefore, when $(m,k)$ is a witness, the shifted pair is a witness exactly
when

$$
 s_p(m,k)+v_p(m+1)-v_p(m+2k)\ge0
 \quad\text{for every prime }p. \tag{5}
$$

Every obstruction in (5) divides $m+2k$, so one factorization gives a complete
shift certificate.

### Odd-$m$ large-prime obstruction lemma

**Lemma.** If $m$ is odd and $(m,k)$ is a witness, every prime obstructing its
natural shift is at most $m$.

**Proof.** Put $x=m+2k$, and suppose $p>m$ divides $x$. The case $m=1$, $p=2$
is impossible because $x$ is odd; hence $p$ is odd. Let $a=v_p(x)\ge1$. For
$q=p^j>m$, Legendre's formula writes the contribution to $s_p(m,k)$ as

$$
D_q=\left\lfloor\frac{x}{q}\right\rfloor
 -2\left\lfloor\frac{x+m}{2q}\right\rfloor. \tag{6}
$$

For $j\le a$, $x/q$ is odd, so $D_q=1$. For $j>a$, write
$x=qt+r$ with $0<r<q$. Both $x$ and $q$ are multiples of $p^a$, hence $r$ is a
nonzero multiple of $p^a$ and

$$
r+m\le q-p^a+m<q.
$$

Equation (6) then gives $D_q=0$ when $t$ is even and $D_q=1$ when $t$ is odd.
Thus $s_p(m,k)=\sum_jD_{p^j}\ge a=v_p(x)$. Also $p\nmid m+1$ in the remaining
odd-$m$ cases. The shifted slack in (5) is therefore nonnegative. $\square$

The parity is load-bearing. When $m$ is even, $x/p^j$ is even for
$j\le v_p(x)$ and those terms contribute zero, allowing the large-prime shift
failures seen in the atlas.

The witness $(m,k)=(3,987)$ is a useful negative control: it is valid, but its
natural shift $(4,986)$ fails exactly at $p=3$. This is consistent with, and
sharp for, the lemma's $p\le m$ conclusion.

## Initial bounded atlas

The exact rectangle $1\le m\le50$, $1\le k\le50{,}000$ contains:

- 12,243 witnesses, all at $m\le10$;
- 96,931 non-witnesses with exactly one obstructing prime;
- 12,243 eligible natural shifts: 6,553 succeed and 5,690 fail;
- among odd $m$, 6,005 of 6,026 shifts succeed; all 21 failures occur at
  $p\le m$;
- among even $m$, 548 of 6,217 shifts succeed; the large-prime mechanism is
  dominant.

### Prime-sensitive atlas repair

The original carry-pattern summaries omitted $p$ from their global keys. That
merged equal-looking carry positions from different bases: 77 old near-witness
groups and 15 old shift groups contained multiple primes.

The dominant old shift group was itself a cross-base collision: count 6,053,
left/right carry positions \([0]\)/\([0]\), first representative \(p=5\), and
1,860 distinct primes in the independently reconstructed split. The regression
suite pins this case explicitly; cross-prime equality is the rule here, not an
edge case. Schema 2 now uses

$$
(p,\text{deficit},\text{left carries},\text{right carries})
$$

and the corresponding prime-first shift key. Exact independent reclassification
changes the number of distinct near patterns from 1,190 to 6,493 and shift
patterns from 337 to 2,249. It identifies 95,330 near events and 6,316 shift
events whose old equivalence class split, for 101,646 exact case records.
`m_rows` and the complete offset screen are byte-identical before and after the
repair. The old and new `reported_instances` values are not comparable: both
report only the top 256 groups, whose granularity changed.

A literal offset screen still leaves 4,104 of 5,687 boundary-comparable failed
shifts unrepaired for every $|\delta|\le16$. This is not a CRT test.

### Structural progress after repair

[`THEOREMS.md`](THEOREMS.md) proves and checks:

- an exact $-1/0/+1$ prime-power zone formula and universal bad-window
  support for every negative level;
- exact prime-power-band CRT classes and uniform local-safe residue classes;
- a Dirichlet construction proving no fixed CRT class can be sufficient;
- the exact decomposition into a small-prime tier and a run of
  compensation-good bad-window integers;
- the full guaranteed odd-prime first-prefix cone, of size
  $\frac{p-1}{2}(\frac{p+1}{2})^{e-1}$;
- a product-displacement lower bound for every repair retaining old blockers;
- a zero-carry CRT system that controls all small primes while evicting any
  prescribed finite set of old large blockers, with its exact class count;
- an adaptive progression that evicts old blockers and visits every uniform
  small-prime-safe class in one exact period;
- unconditional eviction of every attached large prime by translating $k$ by
  $\lceil m/2\rceil$, together with the proof that no smaller uniform offset
  works;
- the powersmooth necessity bound $p^{v_p(w)}<\sqrt{2w}$ at every large prime
  of a compensation-good integer, attained at $w=p^e(p^e+1)/2$;
- the exact $\log 2$ short-cofactor density and the resulting witness density
  ceiling $1-\log 2$ for every fixed $m$;
- the reduction of every natural shift of a witness to one integer: an odd
  source keeps its window and can fail only at $p\le m+1$, while an even source
  adjoins exactly $x/2$;
- sharp short-cofactor and minimal-cofactor prime-power spike thresholds;
- the odd-$m$ CRT-safe shift theorem and the identity equating even-shift
  spikes with compensation of the newly adjoined first bad-window term;
- sharp single-level and smooth-window sufficient criteria.

Finite evidence now separates the mechanisms. All 5,798 CRT-safe odd witnesses
shift successfully. The old fatal-spike lemma explains 4,578/5,669 failed
even shifts; the exact digit formula explains 5,620/5,669, leaving 49 with
only small-prime obstructions. The compensation formula is checked at 50,360
large-prime bad-window events, the short-cofactor corollary at 405,780 factors,
and the band/local CRT formulas at 711 band values and 523,408
unrestricted-high-part extensions.

The band-exact radius-64 cover gives every one of 5,678 boundary-comparable
failed shifts a valid $p\le m+1$ tier candidate. In 5,321 cases the failed
natural target already satisfies that tier. Only 14 first tier candidates are
witnesses; 3,702 are window-blocked before a later witness, and 1,962 have no
witness in the radius. Their first tier windows contain 4,759 single-level
blockers and 1,847 same-prime compensation deficits. Exact small-prime CRT is
therefore almost never the binding constraint in this rectangle.

At the first unlisted index $m=27$, the published $m=26$ witness shifts to a
target killed by $p=5{,}048{,}891{,}644{,}633$. Two one-process segmented
searches check the next 200,000,000 offsets. Single-level primes reject
199,999,988 offsets. Complete checks reject the 12 survivors: the first three
at 14, 6, and 7 large primes; the later nine at 9--17 total primes. Two of the
later cases also fail at $p=7$. The structure certificate checks all 140 large
survivor blockers. Retaining every old blocker forces prime-product
displacements of 30--79 digits, and the earlier zero-carry and adaptive
eviction routes produced 49--98- and 32--81-digit representatives beyond the
retained factorization bound.

That size cost was avoidable. Translating $k$ by exactly $\lceil m/2\rceil$
evicts every attached prime $p>m$ unconditionally, because the translated
window is the next block of $\lceil m/2\rceil$ integers and the two blocks
together span at most $m$ consecutive values. No smaller uniform offset works.
The translated candidate keeps the magnitude of its source, so all 12 survivor
systems give 13-digit candidates that evict all 140 recorded blockers, satisfy
the exact small-prime tier with no search, stay inside the factorization bound,
and are classified completely: each fails only at 9--15 newly entering large
primes. Over $1\le m\le20$, $1\le k\le5{,}000$ the same translation classifies
all 98,706 large-prime-obstructed sources exactly, yielding 1,033 witnesses,
59,641 new-large-prime-only failures, 162 small-prime-only failures, and 37,870
double failures.

The compensation-good criterion removes the residual sieve's proximity
restriction. Two contiguous exact runs exhaust
$50{,}001\le k\le1{,}000{,}050{,}000$. Every one of these $10^9$ candidates
has negative local slack at some $p>27$; the longest observed run of
consecutive $27$-compensation-good integers is 9, short of the required 14.
Together with the atlas, this is an exact finite lower bound on the least
possible $m=27$ witness.

Four new shards apply the same criterion at the shift target's own scale and
extend the certified neighbourhood to the contiguous interval
$5{,}048{,}891{,}644{,}619\le k<5{,}049{,}891{,}644{,}619$, again with no
witness and with longest runs 8, 9, 9, 8. Each shard processed
$2\cdot10^8$ candidates in about 720 s on one reduced-priority process, a
measured $2.8\cdot10^5$ candidates per second. At that rate and the measured
densities, reaching the independence model's predicted $5.1\cdot10^{12}$ would
cost roughly $1.8\cdot10^7$ single-core seconds, so brute force is expensive
but not out of reach; the regression's $2.8\cdot10^{13}$ is about six times
further.

The barrier now has a measured size. Compensation-good density is $0.1237$ at
the $5\cdot10^{12}$ scale against a proved ceiling of $1-\log 2$, the exact
small-prime tier passes for 99.6% of consecutive $k$ there, and observed run
spectra match a term-independence model. The published $m=25$ and $m=26$
witnesses already supply a run of 13 consecutive $27$-compensation-good
integers; the even-to-odd shift adjoins exactly one further term on the left,
and for $m=26$ that term is the prime $5{,}048{,}891{,}644{,}633$, which cannot
compensate.

**Next active proof task:** decide the existence of a length-$\lceil m/2\rceil$
run of $m$-compensation-good integers. Old blockers are no longer part of the
difficulty: translation evicts them at zero cost, the small-prime tier is
almost never binding, and every remaining failure is a newly entering large
prime. The proved powersmooth necessity bound reduces the target to
$\lceil m/2\rceil$ consecutive integers that are simultaneously
$\sqrt{2w}$-powersmooth away from $p\le m$, a correlation problem for large
prime factors of consecutive integers.

## Files

| Path | Role |
|---|---|
| `erdos389.py` | Exact Legendre/Kummer verifier, ratio scanner, floor decomposition, and residue-zone classifier. |
| `atlas.py` | Exhaustive bounded atlas with prime-sensitive pattern identities and offset screen. |
| `compare_atlas_classifications.py` | Independent old/new aggregation re-scan and exact changed-case ledger producer. |
| `shift_repair_analysis.py` | CRT-safe odd shifts, zero-carry moduli, and even prime-power spike certificates. |
| `prime_repair_analysis.py` | Radius-bounded discrimination between original-prime repair and new-prime interference. |
| `zone_analysis.py` | Exhaustive finite certificate for the level-zone, support, and single-level obstruction theorems. |
| `zero_carry_corridor_search.py` | Full or partial small-prime-zero candidate-family search. |
| `residue_control_analysis.py` | Exact band/local CRT classes, fixed-class obstruction examples, bad-window compensation, mirror lifts, and exact shift-spike certificate. |
| `compensation_structure_analysis.py` | Exact local-window separation, short-cofactor spikes, prefix cone, retained-prime displacement, zero-carry blocker eviction, and shift-identity certificate. |
| `adaptive_eviction_analysis.py` | Exact first-hit search along blocker-evicting progressions through uniform local-safe classes, plus a bounded exact new-prime scan. |
| `exact_eviction_certificate.py` | Exhaustive half-window translation eviction, its sharp minimal offset, and exact classification of every translated candidate. |
| `parity_shift_certificate.py` | Exact reduction of every witness shift to its adjoined bad-window term, over a rectangle and all published witnesses. |
| `smooth_density_analysis.py` | Powersmooth necessity checks, two-method short-cofactor counts, measured compensation densities and run spectra, and the labelled growth model. |
| `compensation_run_search.py` | Segmented exact search for runs of compensation-good integers, with full checks of every large-prime-good window. |
| `crt_repair_cover.py` | Radius-bounded separation of exact small-prime tier repair from moving-window blockers. |
| `moving_bad_window_search.py` | Segmented exact single-level sieve plus full checks of every survivor near the first $m=27$ shift target. |
| `bad_window_near_miss_analysis.py` | Full checks and explicit CRT/mirror repairs for closest fatal-window misses and rare sieve survivors. |
| `verify_known.py` | Verifies all 27 OEIS witnesses, 2,000 direct cross-checks, tightness, and two minima. |
| `artifact_io.py` | Same-directory atomic JSON/JSONL artifact writes. |
| `tests/test_erdos389.py` | Arithmetic, key identity, artifact, theorem, and finite-certificate regressions. |
| `THEOREMS.md` | Proofs, corrected failed derivation, finite evidence, and unresolved gaps. |
| `data/oeis_a375071.csv` | Content-addressed local transcription of the OEIS b-file. |
| `data/atlas_m1_50_k50000.json` | Canonical repaired schema-2 atlas. |
| `data/atlas_m1_50_k50000_pre_prime_key_fix.json` | Preserved schema-1 baseline containing the cross-prime merges. |
| `data/atlas_prime_key_classification_diff.json` | Exact split-group summary and extremal classification. |
| `data/atlas_prime_key_changed_cases.jsonl` | All 101,646 reclassified events, SHA-256-pinned by the summary and tests. |
| `data/shift_repair_analysis_m1_50_k50000.json` | Exact bounded CRT/spike analysis. |
| `data/prime_repair_analysis_r64.json` | Exact 5,690-shift prime-aware repair classification. |
| `data/zone_theorem_m1_8_k2000.json` | Exact 16,000-pair zone certificate. |
| `data/zero_carry_corridor_*.json`, `data/partial_zero_carry_*.json` | Exact construction-family positive and negative searches. |
| `data/residue_control_analysis_m1_50_k50000.json` | Exact band, local-safe, compensation, fixed-class, and shift-spike certificate. |
| `data/crt_repair_cover_m1_50_k50000_r64.json` | Exact 5,678-case small-prime-tier versus moving-window repair cover. |
| `data/moving_bad_window_m27_h*.json` | Exact 1M, 20M, and two contiguous 100M moving-window searches near the published natural-shift target. |
| `data/bad_window_near_miss_m27_h*.json` | Exact one-prime and combined-CRT repair analyses for the 1M near misses and 100M sieve survivors. |
| `data/compensation_structure_m1_20_k5000_m27_h*.json` | Exact theorem checks, blocker-eviction experiments, and retained-prime bounds for the 100M and 200M moving certificates. |
| `data/adaptive_eviction_m27_h200000000.json` | Exact first safe progression parameters for all 12 moving-window survivor blocker systems. |
| `data/exact_eviction_m1_20_k5000_m27.json` | Exhaustive eviction certificate, sharpness instances, 98,706 exactly classified translated candidates, and the 12 resolved survivor systems. |
| `data/compensation_run_m27_k5049*_h200000000.json` | Four contiguous exact shards extending the shift-target neighbourhood to $10^9$ candidates. |
| `data/parity_shift_m1_20_k4000_published.json` | Exact shift reduction at 905 rectangle witness sources and all 26 published witnesses. |
| `data/smooth_density_m1_20_k3000.json` | Powersmooth necessity checks, short-cofactor densities against $\log 2$, measured block densities and run spectra, and the labelled growth model. |
| `data/compensation_run_m27_k50001_h100000000.json`, `data/compensation_run_m27_k100050001_h900000000.json` | Two contiguous exact run searches proving the $m=27$ bound through $k=1{,}000{,}050{,}000$. |
| `data/known_witness_verification.json` | Machine-written witness certificates and run summary. |
| `data/atlas_smoke_m1_8_k500*.json` | Repaired and preserved pre-fix smoke artifacts. |
| `run_low_cpu.sh` | One-process, reduced-priority launcher with thread pools pinned to one. |

## Reproduction

From this directory:

```sh
./run_low_cpu.sh python3 verify_known.py
./run_low_cpu.sh python3 atlas.py
./run_low_cpu.sh python3 compare_atlas_classifications.py
./run_low_cpu.sh python3 shift_repair_analysis.py
./run_low_cpu.sh python3 prime_repair_analysis.py
./run_low_cpu.sh python3 zone_analysis.py
./run_low_cpu.sh python3 residue_control_analysis.py
./run_low_cpu.sh python3 crt_repair_cover.py
./run_low_cpu.sh python3 compensation_structure_analysis.py
./run_low_cpu.sh python3 compensation_run_search.py
./run_low_cpu.sh python3 compensation_run_search.py --start-k 100050001 --offset-limit 900000000 --output data/compensation_run_m27_k100050001_h900000000.json
for start in 5049091644619 5049291644619 5049491644619 5049691644619; do
  ./run_low_cpu.sh python3 compensation_run_search.py --start-k $start --offset-limit 200000000 --chunk-size 1000000 --output data/compensation_run_m27_k${start}_h200000000.json
done
./run_low_cpu.sh python3 moving_bad_window_search.py
./run_low_cpu.sh python3 moving_bad_window_search.py --offset-limit 20000000 --chunk-size 1000000 --output data/moving_bad_window_m27_h20000000.json
./run_low_cpu.sh python3 moving_bad_window_search.py --offset-limit 100000000 --chunk-size 1000000 --output data/moving_bad_window_m27_h100000000.json
./run_low_cpu.sh python3 moving_bad_window_search.py --start-k 5048991644619 --offset-limit 100000000 --chunk-size 1000000 --output data/moving_bad_window_m27_h100000000_199999999.json
./run_low_cpu.sh python3 compensation_structure_analysis.py --extended-moving data/moving_bad_window_m27_h100000000_199999999.json --output data/compensation_structure_m1_20_k5000_m27_h200000000.json
./run_low_cpu.sh python3 adaptive_eviction_analysis.py
./run_low_cpu.sh python3 exact_eviction_certificate.py
./run_low_cpu.sh python3 parity_shift_certificate.py
./run_low_cpu.sh python3 smooth_density_analysis.py --short-cofactor-limits 100000 1000000 2000000
./run_low_cpu.sh python3 bad_window_near_miss_analysis.py
./run_low_cpu.sh python3 bad_window_near_miss_analysis.py --offset-limit 100000000 --chunk-size 1000000 --output data/bad_window_near_miss_m27_h100000000.json
./run_low_cpu.sh python3 zero_carry_corridor_search.py
./run_low_cpu.sh python3 zero_carry_corridor_search.py --m-min 8 --m-max 12 --t-limit 10000 --output data/zero_carry_corridor_m8_12_t10000.json
./run_low_cpu.sh python3 zero_carry_corridor_search.py --m-min 27 --m-max 27 --prime-bound 7 --t-limit 20000 --output data/partial_zero_carry_m27_p7_t20000.json
./run_low_cpu.sh python3 zero_carry_corridor_search.py --m-min 27 --m-max 27 --prime-bound 11 --t-limit 20000 --output data/partial_zero_carry_m27_p11_t20000.json
./run_low_cpu.sh python3 -m unittest discover -s tests -v
```

Every artifact records parameters, source/implementation hashes, and its own
resource use. Semantic result fields are deterministic; measured runtimes can
vary. Producers are single-process and atomically replace the final artifact.

## Source status

- [Erdős Problems #389](https://www.erdosproblems.com/389): statement and
  current open-status record, accessed 2026-08-25.
- [OEIS A375071](https://oeis.org/A375071) and its
  [b-file](https://oeis.org/A375071/b375071.txt): listed least witnesses through
  $m=26$, retrieved 2026-08-25.
- Maciej Ulas, *A note on Erdős--Straus and Erdős--Graham divisibility problems
  (with an Appendix by Andrzej Schinzel)*, IJNT 9 (2013), 583--599,
  [DOI 10.1142/S1793042112501497](https://doi.org/10.1142/S1793042112501497).
  The repository metadata and abstract were located; the full article has not
  yet been acquired locally, so detailed literature claims from it remain
  **UNVERIFIED HERE**.
- [FormalConjectures/ErdosProblems/389.lean](https://github.com/google-deepmind/formal-conjectures/blob/main/FormalConjectures/ErdosProblems/389.lean):
  formalized statement and the finite $m=3$, $k=207$ variant.
