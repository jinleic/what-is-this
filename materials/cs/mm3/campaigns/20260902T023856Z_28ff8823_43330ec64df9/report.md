# Campaign report — gate `laderman23-ladder-total` (run 20260902T023856Z_28ff8823_43330ec64df9)

## Question (prereg sha `2af3b859…`, commit `7e3321a`; amendment 1 commit `5e40fc6`)

In the frozen three-stage linear-SLP model, what is the **exact** minimum
addition count for Laderman's own 1976 rank-23 decomposition at its printed
orientation — and how does it compare to the operation count of the scheme as
Laderman printed it?

## Answer — verdict FROZEN-CERTIFIED (`CERTIFIED-EXACT`): the minimum is EXACTLY 62

| stage | certified LB | verified UB | exact |
|---|---|---|---|
| left, `C(U)` | 16 | **16** | yes |
| right, `C(V)` | 16 | **16** | yes |
| output factor, `C(WFac)` | 16 | **16** | yes |
| output stage = `C(WFac) + 14` | 30 | **30** | yes |
| **total** | **62** | **62** | **LB = UB = 62** |

**Laderman's printed basic form recounts to exactly `28 + 28 + 42 = 98`
additions** (recounted twice, independently). The certified minimum for his own
orientation is therefore **36 additions below the printed count**.

**This is not a correction of Laderman.** His paper says explicitly (p.126):
*"Obviously the number of additions in this algorithm could be greatly reduced,
but it is being given in its more basic form."* This campaign **quantifies that
stated slack exactly**: the reducible amount is exactly 36, and 62 is the floor
— no further reduction is possible for this decomposition at this orientation
in this model.

## The reduction and what was actually decided

`total = C(U) + C(V) + C(WFac) + 14` (the `+14 = 23 - 9` transposition gap;
preconditions audited in-run: 23/23 nonzero product rows per side and
`rank(WFac) = 9` over Q).

Nothing was inherited: Phase 1 deliberately claimed no counts. In-run,
`d(U) = d(V) = d(WFac) = 14` with **all three floor verdicts impossible** (2401
DFS states each) ⇒ every side `>= 15` ⇒ `total >= 59`. The decided question was
whether 15 is attainable, i.e. the **aux-1 existential at T = 15**.

**Closure argument (load-bearing).** With `T = 15` and `d = 14`, exactly one
gate creates a non-needed class, so at the moment that auxiliary is created only
inputs and needed classes are available; hence the auxiliary is a signed
pair-sum of `(inputs ∪ tau)`. That universe is finite and was **completely
enumerated and hash-pinned before any feasibility compute**: `|AU| = 408` per
side, sha256 `a1de3df4e5077746` (U), `d3f3894227f939fe` (V),
`13be8045ff2aec8a` (WFac) — re-derived in-run by **two independent enumeration
loop shapes** that agreed with each other and with the pins.

## Decisive evidence

- **Complete census: 3 x 408 = 1,224 dual-instrument instances, ZERO
  admissions.** The memoized subset-DFS and the freshly built slot-availability
  CNF (kissat 4.0.4) agreed on **every single one of the 1,224 instances**;
  index coverage verified to be exactly 0..407 per side with no gaps or
  duplicates. Hence `C(U), C(V), C(WFac) >= 16` and `total >= 62`.
- **All 1,224 infeasibility results carry certificates.** Each UNSAT instance
  produced a kissat DRAT proof converted to LRAT and checked by **both** pinned
  checkers (`drat-trim` `111b0405…`, `lrat-check` `b4bdebfc…`): **1,224/1,224
  accepted with rc 0/0**, recorded per instance in `scan_checkpoint.jsonl`.
  *Retention policy (disclosed):* every certificate was generated and
  dual-checked, and the proof files for the first 5 instances per side plus the
  R3 floor certificate (16 files) are retained in `certificates/`; the rest were
  deleted after passing both checkers, with their `cnf`/`lrat` sha256 and both
  checker return codes retained in the checkpoint. This bounds the frozen dir
  to a sane size while keeping every result reproducible and hash-identified.
- **Upper bound: three explicit 16-gate circuits, each verified independently
  of the search.** The randomized greedy CSE (master seed `20260902`, 240
  restarts/side, jitter schedule fixed in advance) proposed circuits; **240/240
  proposals per side passed exact-expansion verification and 0 were discarded**,
  the best being 16 gates on every side (`witness_U.json`, `witness_V.json`,
  `witness_WFac.json`). Verification recomputes every gate from its operand
  structure over 9 symbolic free inputs and requires all 23 target rows up to
  sign; it shares no code with the search.
- **Output stage constructed, not assumed.** Reverse-mode transposition of the
  verified `WFac` circuit yields an explicit 23→9 stage of exactly **30**
  additions (each accumulation into a non-empty adjoint counted), and it
  reproduces the `WFac` map in **all 207 (23 x 9) entries**.
- **End-to-end.** The assembled scheme (16 + 16 + 30) was expanded symbolically
  over all **81 monomials** `A_i·B_j`: **9/9 outputs equal the 3x3 matrix
  product exactly, zero mismatches** (`assembled_laderman.json`).

## Controls — both directions, all passed

ACCEPT — **A1** frozen triple re-verified (sha `522ba07f4f1784ac…`, Brent
729/729 int + fmpz) and `d`/floor re-derived; **A2 known-true circuit accepted
at its exact addition count**: the frozen `paper55` printed circuits verified by
exact expansion at exactly 13 and 14 gates with the recount `13 + 14 + 28 = 55`;
**A3** frozen landscape row re-verified through the same code path
(`paper55` d/floor = (12,False), (13,False), (13,False)); **A4** 8 pseudorandom
product relabelings per side (seed 17) leave `d` and every floor verdict fixed;
**A5** transposition preconditions; **A6** universe pins via two independent
enumeration loops.

REJECT — **R1** one-addition-deleted plant on a verified 16-gate circuit:
REJECTED (a target class becomes uncovered); **R2** operand-sign-perturbed
plant: REJECTED by exact expansion; **R3** the floor level `T = 14` re-certified
UNSAT through the CNF instrument with a fresh certificate accepted by both
checkers (rc 0/0); **R4** planted wrong `d`-count assertion fired; **R5**
(amended) a creatable chain of 7 classes at `T = 7`, known SAT by construction,
came back feasible from **both** instruments. Cross-solver control: **CaDiCaL
3.0.1** agreed with kissat on the representative instance.

## Independent post-run audit — 8/8 blocks PASS (`postcompute_audit.json`)

Brent re-verified with an independently ordered loop; `d`/floors re-derived;
census completeness re-checked (index coverage, universal agreement, all 1,224
certificate return codes); all three witnesses re-verified by a fresh evaluator;
the output stage **independently re-transposed to 30 additions with all 207 map
entries exact**; an independent end-to-end 81-monomial expansion; 16 retained
certificates replayed through both checkers; verdict arithmetic and the printed
recount (`28 + 28 + 42 = 98`, `98 - 62 = 36`) re-derived.

## Instrument defect (disclosed; caught by the control itself, no claim affected)

**D1 — control R5 was mis-specified in my pre-registration.** The first run
passed every anchor and accept-side control, then aborted at R5 with **no
verdict**. I had pre-registered "a side's classes at `T = 15` + 5 slack slots"
as *trivially schedulable*; that premise is false for this encoding, whose
operand pool is `base = inputs ∪ tau`. Slack slots add no expressive power, so
that instance is UNSAT exactly because the floor is — **the encoder was right
and my control was wrong**. Amendment 1 (committed BEFORE any amended compute)
replaced it with a creatable chain that is known SAT by construction and now
requires agreement from **both** instruments, which is strictly stronger than
the original intent. The aborted verdict and pre-amendment runner are preserved
verbatim (`verdict_preamendment_aborted_R5.json`,
`laderman_ladder_run.py.preamendment`).

## Cost, and the route deliberately not run

Measured **~7.5 CPU-minutes** total (censuses ~78 s/side plus certificate
generation and dual-checking; synthesis ~10 s/side; verification and assembly
seconds) against a pre-registered **2.0 CPU-hour** cap. No budget change, no
sampling: every claim rests on a complete census.

**Not run, quantified:** the aux-2 census at `T = 16` (~90,000 pairs/side,
measured ~2.0 CPU-h/side, ~6 CPU-h total). It was **never needed** — the aux-1
census already gave `LB = 16` per side and the verified 16-gate circuits met it
exactly, so the ladder closed one level earlier than the pre-registered budget
allowed for. It is recorded here as a quantified route that was not taken and
from which nothing is claimed.

## Frontier statement (bounded, precise)

In the three-stage linear SLP model over `{-1,0,1}` (inputs free; gate `= x ± y`
costs 1; sign changes and copies free; output stage counted with the
transposition bound), for the factor triple transcribed verbatim from
Laderman's 1976 paper at its printed (sigma^0, all-monomial) orientation:
**the minimum addition count is exactly 62 = 16 + 16 + 30**, attained by an
explicit circuit whose 9 outputs are verified equal to the 3x3 matrix-product
bilinear forms by exact symbolic expansion over all 81 monomials. The printed
basic form uses 98; the reducible slack the paper itself advertises is **exactly
36 additions, no more**.

This says nothing about other orientations of Laderman's tensor, other
decompositions, non-ternary alphabets, `GL(3,Q)`/`GL(3,Z)` sandwiches, other
actions or counting conventions, or the global sub-55 question. Laderman's 62 is
**above** the record 55 (`paper55`), as expected for a different decomposition —
it is not a competitor for the record, and no claim about the record is made.

## Exact fixed-orientation ladder after sessions 9-15

| decomposition / orientation | published | certified exact minimum |
|---|---|---|
| `paper55` sigma^0 | 55 | **55** (record) |
| `paper55` sigma^1 | (57 best known) | **55** |
| `paper55` sigma^2 | (57 best known) | **55** |
| `sun56` sigma^0 | 56 | **56** (published optimal) |
| `mws59` sigma^0 | 59 | **58** (published +1) |
| `stapleton60` sigma^0 | 60 | **60** (published optimal) |
| **`laderman23` sigma^0** | **98** (explicitly unoptimized) | **62** |

## exactly-one verdict

**FROZEN-CERTIFIED** — `C(U) = C(V) = C(WFac) = 16` exactly, output stage 30
exactly, **total exactly 62 (LB = UB)** for Laderman's fixed orientation; lower
bound from complete floor DFS plus a complete hash-pinned 1,224-instance
dual-instrument census with 1,224 dual-checked certificates, upper bound from
explicit circuits verified independently of the search and end to end over all
81 monomials; Laderman's printed 98 exceeds the true optimum by exactly 36.
