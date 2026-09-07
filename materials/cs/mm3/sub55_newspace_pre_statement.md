# Pre-statement — `sub55-newspace`: the full ternary orientation orbit of `laderman23`, the last unswept named slice of the sub-55 question

**Created 2026-09-03, agent `Mm3Sub55`, BEFORE any search compute of this
campaign.** Committed path-scoped before `campaign.py init`; the minted run
dir receives a byte-identical copy with source commit + sha256 recorded.

## 0. State inherited from frozen evidence (read, not recomputed)

- The fixed-orientation ladder is settled and NOT re-litigated here:
  `paper55` sigma^0/1/2 = 55, `sun56` = 56, `mws59` = 58, `stapleton60` = 60,
  `laderman23` = 62 (campaigns of sessions 5-15; latest
  `20260902T023856Z_28ff8823_43330ec64df9`, FROZEN-CERTIFIED).
- Sessions 5-8 swept the FULL ternary-unimodular orientation orbit
  $(P,Q,R)\in T^3$ (|T| = 6960) x sigma-orbit {0,1,2} for exactly FOUR factor
  triples: `paper55`/`perminov58`, `sun56` (session 7), `mws59`,
  `stapleton60` (session 8). Zero <=54 anywhere; all certified >= 55.
- `laderman23`'s verified factor triple did not exist until session 14
  (`20260902T022208Z_8c279ad1_d0276dfee63d`, factors Brent-verified 729/729
  over Z in int and fmpz, triple sha256 `522ba07f4f1784ac...`). Session 15
  decided ONLY its printed sigma^0 orientation (exactly 62) and its scope
  sentence states it "says nothing about other orientations of Laderman's
  tensor". No frozen artifact contains any d-count, floor verdict, or bound
  for any sigma^1/sigma^2 orientation or any non-monomial sandwich of
  `laderman23`. Session 8 named exactly this follow-up: "lock a primary ...
  then pre-register its ternary-unimodular off-diagonal census if and only if
  the factor alphabet fits this model" — the alphabet fits (all-ternary
  factors, verified), and this is that pre-registration.
- Frozen reductions this campaign reuses (all PROVED, not assumed):
  - **Monomial-transfer theorem** (session 5): $C(F)=C(gF)$ for every signed
    permutation $g$ of the 9 input coordinates; the 48^3 monomial sandwich
    core is cost-invariant.
  - **145-node data factorization** (sessions 7-8): $T=$ 6960 = 145 data
    nodes x 48 right-monomials; every sandwich predicate depends only on the
    data triple, so each data triple represents exactly 48^3 = 110,592
    sandwich triples.
  - **Automorphism family** (session 7, proven and controlled): in standard
    matrix semantics, blockwise per 3x3 block,
    $U'=P^{-1}UQ^{-\mathsf T}$, $V'=Q^{\mathsf T}VR^{-\mathsf T}$,
    $W'=P^{\mathsf T}WR$ preserves the tensor (Brent 729/729) for every
    invertible $P,Q,R$; the all-ternary predicate and all bound quantities
    are decided on the 145^3 data triples.
  - **Summand-order invariance** (session 5): $d(F)$ is a cardinality of a
    set of sign classes; permuting/flipping the 23 summands cannot change it.
  - **Transposition accounting** (frozen model of sessions 5-15):
    total $= C(U)+C(V)+C(W_{\text{out}})$, lower-bounded by
    $\sum_{\text{sides}}(d+[{\rm floor\ impossible}]) + 14$ where
    $+14=23-9$ is the output stage's transposition gap, valid when all 23
    products are used and all 9 outputs are nonzero (activity precondition,
    audited per instance in this campaign — a per-instance strengthening of
    sessions 7-8, which audited at sigma-class level).

## 1. The three candidate spaces, weighed (as assigned)

- **(a) Full orientation orbit of a fixed decomposition** — chosen, applied
  to `laderman23`. This is the unique remaining *named-public* slice: every
  other public rank-23 factor triple already has its full orbit frozen
  (sessions 5-8). New relative to all frozen work: no frozen artifact
  decided any orientation of `laderman23` other than sigma^0.
- **(b) Data-triple variation at fixed factors** — subsumed. At fixed
  factors with only the monomial sandwich, the data-triple space IS the
  sigma-orbit (3 points); for `laderman23` those 3 points are inside the (a)
  census. For the other four decompositions the general data-triple census
  (all survivor triples, 0 at <=54) is already frozen (sessions 7-8: 5,700 +
  4,510 + 5,796 + 4,800 decided). A (b)-only campaign would decide strictly
  less than (a).
- **(c) Decomposition-independent stagewise lower bound** — rejected with a
  closed-form impossibility argument. The sharp universal per-side bound is
  $C(F)\ge \operatorname{rank}(F)-9$ (span/counting argument), and it is
  TIGHT: for any $r\le 23$ there are 23-row factor matrices of rank $r$
  computable in exactly $r-9$ additions (the 9 unit directions free; each
  further row $e_1+e_i$ one addition). Hence for EVERY rank-23 decomposition
  with factor-rank profile $(r_u,r_v,\cdot)$,
  $$\text{total}\ \ge\ (r_u-9)+(r_v-9)+14 = r_u+r_v-4\ \le\ 23+23-4=42,$$
  strictly below 54. No rank-profile stagewise argument can ever certify a
  global >=55 bound; reaching 55 would require a genuinely new global
  mechanism (a publishable-scale result), out of scope for one campaign.
  Route (c) is therefore frozen here as quantifiably non-decisive, with this
  arithmetic as its measured cost; it MUST NOT be re-attempted at this gate
  without a new mechanism named in a fresh prereg.

## 2. Fixed question and falsifiable outcome

**Question.** In the frozen three-stage linear-SLP model (alphabet
{-1,0,1} coefficient factors; additions = subtractions, counted; copies and
sign changes free; total $=C(U)+C(V)+C(W_{\rm out})$), does ANY orientation
of the `laderman23` factor triple within the proven automorphism family
$(P,Q,R)\in T^3$, $|T|=6960$ ternary unimodular matrices, composed with the
sigma-orbit {Id, sigma, sigma^2}, $\sigma:(U,V,W)\mapsto(V,W^{\mathsf T},
U^{\mathsf T})$, admit a certified three-stage total lower bound <= 54?

**Space, exactly counted (fixed now).** Raw scheme-instances:
$6960^3\times3 = 1{,}011{,}460{,}608{,}000$. Enumerated through the proven
factorization as $145^3\times3 = 9{,}145{,}875$ data-level points; every
all-ternary data triple (a "survivor") is decided exactly. Expected survivor
count 4,000-7,000 by analogy with the four swept decompositions (4,510-5,796)
— **scouting estimate only, never a result**.

**Falsifiable outcome.** A single survivor with certified total LB <= 54
falsifies the negative branch. Note an LB <= 54 is a *candidate*, not a
scheme; see verdict rule.

## 3. Closure and symmetry reductions (PROVED — what collapses the space)

1. Right-monomial fiber: 6960 = 145 x 48; decisions at data level; factor
   48^3 per data triple (session 7, proven + direct/table controls re-run
   here on `laderman23`).
2. Monomial-transfer theorem: within the monomial core, cost-invariant; the
   all-monomial data triple of each sigma class is the class's representative
   anchor (session 5, frozen; its corollary is re-machinery-checked here).
3. Summand relabeling invariance of $d$ (session 5; re-checked here with 8
   pseudorandom relabelings per tested orientation).
4. Automorphism family preserves Brent validity (session 7's symbolic
   delta-collapse + battery; here the family is re-anchored on `laderman23`:
   every sigma class of the base triple must be Brent 729/729 under the
   family at monomials, and the planted wrong-action variants must fail).
5. sigma-orbit: NO proven reduction collapses it; all three classes are
   enumerated and decided exhaustively. (sigma is a data-level cyclic
   automorphism not expressible as a block-diagonal sandwich of the same
   base triple; the three slices are decided independently.)

**What remains after reductions (quantified):** $1.011\times10^{12}$ raw
instances collapse to $9{,}145{,}875$ data-level points; the all-ternary
predicate admits the survivors (in-run count, hash-pinned pair tables);
every survivor is decided by exhaustive exact instruments. Nothing is
sampled.

## 4. Instruments and versions (pinned before use)

- Interpreter: Python 3.14.3 in dedicated venv `~/.venvs/mm3sub55`
  (python-flint 0.9.0, numpy 2.4.3 — byte-identical version set to the
  frozen session-15 toolchain; the pre-existing environment had lost
  python-flint, disclosed here; the venv restores the pinned versions).
- Frozen repository instruments, imported unmodified (file sha256s recorded
  in-run into `versions.txt`): `src/new_decomp_offdiag.py` (library:
  pair tables, survivors, side decisions, controls), `src/gate_b_floor.py`
  (frozen DFS, cross-instrument), `src/gatec_sweep.py` (frozen d-counter),
  `src/gatec_decomps.py` (control loaders), `src/verify_anchors.py`,
  `src/tensor_data.py`, `src/gate_a.py` (Brent batteries).
- SAT layer (NEW relative to sessions 7/8, which were DFS-only):
  kissat 4.0.4 (`/opt/homebrew/bin/kissat`) on the slot-availability CNF at
  $T=d$ (session-15 encoding `build_ext_cnf`, dimension N=9, base universe =
  the 9 input directions plus the needed classes; both-direction slot-0
  units); certificates DRAT -> LRAT via drat-trim, accepted by BOTH pinned
  checkers drat-trim and lrat-check (frozen 2e3b2dc snapshots, sha256
  `111b0405566d55629f5d...` / `b4bdebfcc40da664...`, copied byte-identical
  from the session-15 run dir `tools/`). Cross-solver CaDiCaL 3.0.1 on
  deterministic representatives.
- Process: `nice -n 10`; single-thread BLAS/OMP env pinned; RLIMIT_CPU soft
  raised to hard at the top of every long stage; PYTHONDONTWRITEBYTECODE=1
  and `sys.dont_write_bytecode = True` before importing anything from frozen
  campaign directories.

## 5. Hash-pinned universes (pinned NOW, before any compute)

1. `laderman23` factor triple:
   `campaigns/20260902T022208Z_8c279ad1_d0276dfee63d/laderman23_factors.json`,
   file sha256 pinned in-run; in-run recomputation of the canonical triple
   hash (json sort_keys of {U,V,W}) MUST equal the frozen
   `522ba07f4f1784ac8c31631030027acfe4ee72d166f315ac15044520e7b81b84`.
   Frozen file read read-only.
2. Data-node universe: the 145 representatives, enumerated by TWO independent
   loop shapes (frozen `build_data_orbits` and a fresh brute force over all
   $3^9=19{,}683$ ternary matrices -> 6960 unimodular -> right-monomial
   orbits); both canonical serializations hashed; hashes MUST be equal and
   are pinned in the report.
3. Pair tables (3 sides x 3 sigma, 145x145 boolean each), survivor list,
   and per-pair decision tables: sha256-pinned as written.
4. Known-true circuit witnesses (ACCEPT controls): `paper55` printed SLPs
   (`tensor_data.LEFT_SLP/RIGHT_SLP` + output network) at exactly 13+14+28;
   session-15 `witness_U.json`/`witness_V.json`/`witness_WFac.json` at
   exactly 16/16/30. File hashes pinned in-run.
5. Startup anchor table = frozen `EXPECTED_CLASSES` (12 rows, four control
   decompositions x 3 sigma) plus the `laderman23` sigma^0 row (14,14,14),
   floors impossible, total LB 59, from the frozen session-15 report.
6. Certificate ledger: per instance, sha256 of CNF and LRAT proof, return
   codes of drat-trim conversion and both checker verifications, solver
   seconds.

## 6. Certificate policy (the load-bearing upgrade over sessions 7/8)

Every floor-impossible side instance among decided pairs is decided by BOTH
instruments: (i) complete memoized subset-DFS (two independent
implementations: this campaign's own and frozen `gate_b_floor`), and (ii)
kissat on the slot-availability CNF at $T=d$. The verdicts MUST agree on
every instance. Every infeasibility carries a DRAT proof converted to LRAT
and accepted by BOTH pinned checkers (rc 0 both), recorded in the ledger.
Retention: ALL `.lrat` proofs are retained in the run dir; CNFs are retained
for 12 deterministic representatives per side (lexicographically first per
sigma class) plus all anchor/control instances, and are otherwise exactly
regenerable from the ledger's per-instance target rows + the frozen builder
(builder hash pinned). CaDiCaL must agree with kissat on all representatives.
kissat per-instance timeout 600 s, one retry at 1800 s. If any
floor-impossible instance lacks a dual verdict or accepted certificate after
retry, the negative branch is UNAVAILABLE and the verdict is
FROZEN-INCONCLUSIVE naming exactly the uncovered instances.

## 7. Budget (fixed now; measured cost reported either way)

| phase | cap |
|---|---|
| A: anchors + controls | 0.5 CPU-h |
| B: data-node universe + pair tables (9 tables) + survivor enumeration + table controls | 1 CPU-h |
| C: DFS decisions (all positive pairs, 3 sides x 3 sigma) | 4 CPU-h |
| D: CNF dual layer + certificates | 12 CPU-h |
| total wall | 24 h |

Observed frozen rates (session 7: 929 s + 809 s for two decompositions'
complete DFS censuses; session 15: 1,224 CNF instances + certificates in
~7.5 CPU-min) make full coverage likely inside 2-4 h. If any cap binds: NO
partial negative — verdict FROZEN-INCONCLUSIVE with the exact decided prefix
named and the undecided remainder explicitly out of scope. A stage that dies
(SIGXCPU exit 152, OOM, or stall) is relaunched only if its artifact state
allows exact resume; otherwise the same FROZEN-INCONCLUSIVE rule applies.
Polling: the runner polls its own progress artifacts; no blocking wait
longer than a stage's expected wall time.

## 8. Controls — both directions, in-run, AFTER init, BEFORE any claim

**ACCEPT (must pass; abort otherwise):**
- A1: `laderman23` base triple Brent 729/729 in BOTH arithmetics (int,
  fmpz) after byte/hash verification of the frozen factor file.
- A2: startup anchors — all 12 frozen EXPECTED_CLASSES rows reproduce
  (d, floor, total 55/55/56/58) through both the campaign DFS and frozen
  `gate_b_floor`; `laderman23` sigma^0 reproduces (14,14,14)/impossible/59
  through both instruments.
- A3: known-true circuits accepted at their exact addition counts — the
  verifier accepts `paper55`'s three circuits at 13/14/28 (recount 55) and
  the session-15 witnesses at 16/16/30 (recount 62), each by exact symbolic
  expansion over its full target set (independent of any search).
- A4: relabeling invariance — 8 pseudorandom relabelings of the 23 products
  leave (d, floor, LB) fixed on `laderman23` sigma^0 and on 2 survivor
  triples per sigma class.
- A5: universe pins — the two independent 145-node enumerations hash-equal;
  48^3 touchstone (110,592) and 145^3 counts asserted; direct-vs-table
  full-row subcube (145^2) + 200-random concordance + right-monomial fiber
  predicate checks on `laderman23`.
- A6: activity audit — per sigma class and per survivor triple: all 23
  products used, all 9 outputs nonzero (precondition for the +14 gap);
  any inactive triple drops the gap (conservative LB) and is reported.

**REJECT (must fail as planted; abort otherwise):**
- R1: gate-deleted plant — delete one gate from a known-true `paper55`
  circuit; the exact-expansion verifier must REJECT.
- R2: operand-perturbed plant — flip one operand's sign in a known-true
  circuit; must REJECT.
- R3: planted floor instances — `laderman23` U-side sigma^0 at $T=14$ must
  be UNSAT in BOTH instruments with a fresh dual-checked certificate; a
  known-creatable chain instance (session-15 amended-R5 pattern, e.g. class
  $e_1+e_2$ at $T=1$ from inputs) must be SAT in the CNF instrument.
- R4: planted-count anchors — the all-monomial data triples must total
  exactly 55 (paper55 control), 56 (sun56), 56 (mws59), 58 (stapleton60),
  59 (laderman23 LB anchor); any mismatch aborts before histograms are
  written.
- R5: wrong-action plants (frozen session-7 pattern) — the transpose-as-
  inverse and wrong-$W'$ sandwich variants must FAIL the ternary/Brent
  predicates exactly as frozen.

No transposition-derived witness is used anywhere: the only witnesses in
this campaign are the frozen session-15 circuits, re-verified by direct
exact expansion, and DFS floor schedules, replayed gate by gate.

## 9. Verdict rule (fixed before compute) — exactly one terminal verdict

- **FROZEN-NEGATIVE** (expected): every survivor data triple of all three
  sigma classes decided, 100% dual-instrument agreement, 100% certificate
  coverage (section 6), and minimum certified total LB >= 55 (zero rows
  <= 54). The claim is EXACTLY: *no ternary orientation of the `laderman23`
  factor triple within the proven automorphism family $(P,Q,R)\in T^3$ and
  the sigma-orbit has a certified three-stage total lower bound <= 54; the
  certified LB histogram and minimum over the survivor set are as reported;
  combined with frozen sessions 5-8, the complete named landscape of five
  public rank-23 factor triples (paper55/perminov58, sun56, mws59,
  stapleton60, laderman23) is now certified >= 55 over its full ternary
  orientation orbit, with zero <=54-certifiable points.* An LB is a lower
  bound only: no synthesis, no circuit, no "55 is optimal" claim follows.
- **LIVE / escalation branch**: any survivor with certified LB <= 54 halts
  the runner BEFORE recording the offending row (RecordFlag, frozen
  pattern); agent escalates to Main via hub BEFORE any claim is written; a
  sub-54 LB is a candidate for a separately pre-registered synthesis
  campaign, never a record.
- **FROZEN-INCONCLUSIVE**: any budget cap binds, any instrument disagreement
  survives retry, or certificate coverage is incomplete — with the exact
  decided set, measured cost, and the precise uncovered remainder named.
- Control failure: abort with NO verdict (preserve outputs verbatim),
  amendment with both prereg hashes before any amended rerun.

## 10. What this campaign may NEVER claim

No claim about: any decomposition outside the five named public ones
(Smirnov's other schemes, Schwartz–Vaknin, unpublished, F_2-only — F_2
validity does not imply Brent over Z); non-ternary alphabets;
GL(3,Q)/GL(3,Z) sandwiches beyond ternary unimodular (P,Q,R); actions
outside the proven automorphism family; anti-cyclic (U,V)-swaps computing
BA; any upper bound beyond the frozen witnesses; the global question of
whether ANY rank-23 scheme of <= 54 exists — this campaign completes the
named slice and nothing more. A budget-bound negative, if it were to occur
on a pre-registered sub-slug (it may not, per section 7's no-partial rule),
would name its exact space and measured cost and nothing else.

## 11. Reproduction command (fixed)

```
~/.venvs/mm3sub55/bin/python <run_dir>/sub55_newspace_run.py --run <run_dir>
```

with `PYTHONDONTWRITEBYTECODE=1`, from `cs/mm3`. The runner aborts on any
control failure, any instrument disagreement, and any <= 54 row (before
recording it).
