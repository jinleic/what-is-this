# Pre-statement — `sixth-source-or-orbit-min`

**Created 2026-09-04, agent `Mm3SixthSource`, BEFORE any claim-relevant
compute of this campaign.** Committed path-scoped (cs/mm3 only) before
`scripts/campaign.py init`; the minted run dir receives a byte-identical
copy with source commit + sha256 recorded (PREREG_PROVENANCE.json).

This campaign has TWO preregistered branches, exactly one of which is
executed and closed; the branch is selected by the phase-B source audit
below, mechanically, per the rules fixed here. Exactly ONE terminal verdict
is recorded (section 9).

## 0. State inherited from frozen evidence (read, not recomputed)

- Fixed-orientation ladder, all exact (sessions 9-15):
  `paper55` sigma^0/1/2 = 55, `sun56` = 56, `mws59` = 58, `stapleton60` = 60,
  `laderman23` = 62 (its printed basic form recounts to 98).
- Session 16 (run `20260904T015105Z_49ed8737_988c6df478f1`,
  FROZEN-NEGATIVE) decided the full ternary orientation orbit of
  `laderman23`: 34,944 admitted data triples; per sigma class the
  certified-LB histogram has minimum 59 (unique, the all-monomial data
  triple (125,125,125), d=(14,14,14), floors impossible) and exactly 12
  rows at 61 (all floors impossible; per row two sides at d=15 and one at
  d=14), zero rows at 60, zero rows <= 54. Frozen decision artifact:
  `campaigns/20260904T015105Z_49ed8737_988c6df478f1/decisions_summary.json`.
- The named landscape (five public rank-23 factor triples x full ternary
  orientation orbit x sigma) is closed at the LB layer. The two
  logically-next items, both named in the session-16 report, are exactly
  this campaign's two branches:
  (i) a sixth public rank-23 decomposition, source-locked and decided on
  its own registered orientation slice;
  (ii) `laderman23-orbit-min` — the exact minimum of the `laderman23`
  ternary orientation orbit, ∈ {59, 61, 62}.
- Machinery reused (all frozen, hash-checked in-run before use, never
  modified): `src/new_decomp_offdiag.py` (pair targets, data orbits, Brent
  batteries, EXPECTED_CLASSES), `src/gate_b_floor.py` (prep / subset_dfs /
  canon — the frozen DFS), `src/gatec_sweep.py` (d-count, sigma orbit),
  `src/tensor_data.py`, `src/gate_a.py`, `src/verify_anchors.py`; the
  session-15 aux-1 instruments (`au_universe`, `build_ext_cnf` CNF +
  kissat + DRAT->LRAT dual-checker certificates, reverse-mode
  transposition, exact-expansion circuit verification) copied verbatim
  from `campaigns/20260902T023856Z_28ff8823_43330ec64df9/
  laderman_ladder_run.py` (file sha256 recorded in-run); the checkers
  copied byte-identical from that run dir's `tools/` (drat-trim
  `111b0405566d...`, lrat-check `b4bdebfcc40da664...`, the frozen
  2e3b2dc snapshots); kissat 4.0.4 and CaDiCaL 3.0.1 from
  /opt/homebrew/bin (versions.txt). Python 3.14.3 venv `~/.venvs/mm3sub55`
  (python-flint 0.9.0, numpy 2.4.3).
- Frozen hashes this campaign binds to (verified in-run, phase A):
  `laderman23_factors.json` triple sha256
  `522ba07f4f1784ac8c31631030027acfe4ee72d166f315ac15044520e7b81b84`;
  145-node universe serialization hash prefix `bed89ca17f5bd868`;
  session-16 decision artifact sha256 pinned in-run at first read.
- Process: `nice -n 10`; single-thread BLAS/OMP; RLIMIT_CPU soft raised to
  hard at the top of every long stage; PYTHONDONTWRITEBYTECODE=1 and
  `sys.dont_write_bytecode=True` before any import from frozen campaign
  directories. No writes into any frozen campaign directory (frozen files
  are opened read-only; copied instruments live in THIS run dir).

## 1. Branch selection — the phase-B source audit (bounded, rules fixed now)

**Named candidate sources, in priority order (the audit is restricted to
exactly these; discovering a NEW source mid-run is out of scope and would
need a fresh prereg):**

- **S1 Smirnov 2013** — A. V. Smirnov, "The bilinear complexity and
  practical algorithms for matrix multiplication", Comput. Math. Math.
  Phys. 53(12):1781-1795 (2013); Russian original Zh. Vychisl. Mat. Mat.
  Fiz. 53(12):1970-1984; DOI 10.1134/S0965542513120129. The only public
  rank-23 3x3 scheme from outside the 2025-26 record lineage with
  "fewest known nonzeros" (84 additions naive).
  - Data path A (primary): mathnet.ru full text of the paper
    (paperid zvmmf9955), fetched and hash-pinned; scheme transcribed from
    the printed tables.
  - Data path B (secondary digitization, README-cited by its repository
    to the primary): `github.com/arbenson/fast-matmul` files
    `codegen/algorithms/smirnov333-23-139` and
    `codegen/algorithms/smirnov333-23-128` (27x23 ternary coefficient
    arrays, `#`-separated, layout = 9 U-rows, 9 V-rows, 9 W-rows x 23
    products), fetched at a pinned commit and hash-pinned.
  - Lock rule: a path locks S1 iff the fetched data, transcribed into the
    frozen 23x9 U/V/W convention, passes Brent 729/729 over Z in BOTH
    arithmetics (int, fmpz). Path A and path B are locked (or fail)
    independently; either suffices for a candidate, with the provenance
    level recorded honestly (path A = CITED-DEPENDENCY primary;
    path B = CITED-DEPENDENCY secondary digitization).
- **S2 Schwartz-Vaknin 2023** — O. Schwartz, N. Vaknin, "Pebbling game and
  alternative basis for high performance matrix multiplication", SIAM J.
  Sci. Comput. 45(6):C277-C303 (2023), DOI 10.1137/22M1502719; author
  public PDF (`cs.huji.ac.il/~odedsc/papers/SISC23-Pebbling.pdf`),
  fetched and hash-pinned. Lock rule: iff the PRIMARY TEXT (or its pinned
  official artifact, bounded artifact search <= 30 min) contains a
  complete rank-23 3x3 factor triple (all 23 products' coefficients).
  Recon finding recorded BEFORE this prereg: the paper's bilinear phase
  content is 7x7 alternative-basis Strassen; no 3x3 rank-23 factor data
  was seen in the text. The in-run audit re-checks this first-hand and
  records the outcome; alternative-basis schemes are NOT convertible
  candidates here (their data is not printed; constructing a scheme
  myself is not source-locking).
- **S3 Heule-Kauers-Seidl 2019** — M. Heule, M. Kauers, M. Seidl, "New
  ways to multiply 3x3-matrices", arXiv:1905.10192; J. Symbolic Comput.
  104:899-916 (2021). Lock rule: iff a complete factor triple of one of
  their new rank-23 schemes is retrievable from the paper's pinned
  ancillary/official data within <= 30 min of bounded search. (Their
  schemes are the known family of >13,000 mutually inequivalent rank-23
  schemes; any ONE locked triple suffices for the branch-A question.)
- **S4 Benson-Ballard framework schemes** — A. R. Benson, G. Ballard,
  "A framework for practical parallel fast matrix multiplication",
  arXiv:1409.2908 / PPoPP 2016, repository `github.com/arbenson/
  fast-matmul` files `codegen/algorithms/fast333-23-125`,
  `fast333-23-152`, `fast333-23-221` (one or more; each independently a
  candidate), fetched at a pinned commit and hash-pinned.
- **S5 (control-only) the five frozen triples** — paper55, sun56, mws59,
  stapleton60, laderman23 (frozen factor artifacts), which MUST pass
  through the same ingest/verify path as the candidates (ACCEPT
  controls), and against which the equivalence test is calibrated.

**Audit budget: 2.0 h wall total for all fetches, transcriptions,
verifications, and the equivalence battery. If the budget binds, the
un-locked remainder is recorded as not-locked and the branch rule below
applies with exactly the candidates locked so far.**

**Branch rule (mechanical):** if some candidate locks S1/S3/S4 (i.e.,
Brent 729/729 over Z in both arithmetics) AND passes the
non-equivalence test of section 2 against ALL FIVE frozen triples ->
**Branch A** (sixth-source orientation search, section 4). Otherwise
(including: S2-only, all candidates equivalent to a frozen triple, all
fetches failing, budget binding) -> **Branch B** (`laderman23-orbit-min`,
section 5). The audit's full negative evidence (per source: what was
fetched, its sha256, why it failed to lock) is frozen either way; a
negative branch is NEVER recorded without this evidence file.

## 2. Equivalence test (exact, fixed now; decides "not merely a
signed/permuted/sigma copy")

For integer triples in the frozen 23x9 convention, define the free action
group G = (sigma^j, j in {0,1,2}) x (g, any of the 48 signed 3x3
permutation matrices) x (h, any of the 48), acting by
(U,V,W) |-> sigma^j( gU-blocks, gV-blocks, hW-blocks ) with the 9 input
coordinates of U and V conjugated by g and the 9 output coordinates of W
by h, composed with per-product rational normalization
(U_k, V_k, W_k) |-> (U_k/u0, V_k, u0 W_k) where u0 is U_k's first
nonzero coefficient in the fixed 9-coordinate scan order, and with an
arbitrary permutation of the 23 products. Canonical form = the sorted
list of the 23 normalized product triples, serialized. Candidate C is
EQUIVALENT to frozen F iff canon(C, x) == canon(F, e) for some x in G
(e = identity), decided by exhaustive evaluation over G (3 x 48 x 48 =
6,912 canonicalizations per triple, exact integer arithmetic).
Distinctness from ALL FIVE frozen triples is required for Branch A.
Rationale recorded: G contains per-column signs, product permutations,
sigma, the T involution (coordinate permutation [0,3,6,1,4,7,2,5,8]),
entry transposition (as a sigma image), and monomial basis relabelings —
strictly more than the session-2 Sun/paper55 test used.

## 3. Controls, both directions, in-run, after init, before any claim

**ACCEPT (all must pass; any failure aborts with NO verdict):**
- A1: frozen `laderman23` triple re-verified 729/729 Brent over Z in int
  AND fmpz through this campaign's ingest path.
- A2: the five frozen triples re-verified 729/729 (int) and their
  EXPECTED_CLASSES d-counts ((12,13,13)... etc., 12 rows) reproduced
  through the frozen d-counter on the sigma^0 rows.
- A3: session-15 witnesses (witness_U/V/WFac.json, 16/16/16) accepted by
  exact expansion at exactly 16 gates against the frozen sigma^0 rows;
  paper55 printed SLP accepted at 13/14 (recount 55).
- A4: equivalence-test calibration — a planted scrambled copy of sun56
  (random product permutation x random g x per-product sign flips) MUST
  be judged EQUIVALENT to sun56; sun56 vs paper55 and laderman23 vs
  paper55 MUST be judged NOT equivalent (frozen session-2 facts).
- A5: universe pins — the 145-node representatives rebuilt two
  independent ways (frozen `build_data_orbits`; fresh brute force over
  3^9), serializations hash-equal to each other and to the frozen
  `bed89ca17f5bd868` prefix; session-16 decision artifact hash-pinned at
  first read.
- A6: instrument preflight — one known-UNSAT floor instance
  (laderman23 sigma^0 U at T=14, frozen) returns UNSAT in the DFS and the
  CNF instrument with a fresh dual-checked certificate; one known-SAT
  creatable-chain instance (session-15 amended-R5 pattern) returns SAT in
  the CNF instrument; certificate checkers replay the frozen session-15
  retained proof `WFac_aux4.lrat` with rc 0 (both).

**REJECT (must fail as planted; abort otherwise):**
- R1: sign-flip one entry of a candidate's fetched data -> Brent must
  break (>= 1 identity).
- R2: delete one product row -> Brent must break.
- R3: planted-count anchors — the three sigma^0 monomial LB anchors
  55/56/58 (paper55/sun56/mws59) and laderman23 sigma^0 LB 59 must
  reproduce through this campaign's own decision arithmetic on the
  monomial triples; mismatch aborts before any census row is trusted.
- R4: equivalence plants — the scrambled sun56 copy (A4) judged
  equivalent, and a corrupted scrambled copy (one product triple altered)
  judged NOT equivalent.

## 4. Branch A — sixth-source orientation search (registered space)

**Space (fixed now):** the newly locked factor triple D6, under the proven
automorphism family (P,Q,R) in T^3 (|T| = 6960 ternary unimodular) x
sigma-orbit {0,1,2} — exactly the session-16 space with D6 in place of
`laderman23` — IFF D6 is all-ternary (verified in-run). If D6 is integer
but not ternary, the registered space degrades to the sigma-orbit at the
monomial triple only (3 orientations), with per-side exact aux-1 censuses
as in section 5's instrument set. Ternary/non-ternary is decided by the
data, before the census starts, and recorded.

**Instruments:** frozen session-16 protocol: d-count + floor verdict per
side (frozen DFS + gate_b_floor cross-check on deterministic samples),
slot-availability CNF at T=d (kissat) dual on every positive pair-side
instance, DRAT->LRAT dual-checker certificates on every floor-impossible
instance, 145-node factorization with two-shape universe pin, direct-vs-
table controls, activity audit per triple, relabeling invariance 8x.

**Verdict content:** certified-LB histogram and minimum over all admitted
data triples of D6's full registered space; escalation: any triple with
certified LB <= 54 halts the runner BEFORE recording the row and is
escalated to Main via hub before any claim (a sub-54 LB is a candidate,
never a record). Budget: 14 CPU-h compute + 2 h wall audit phase
already spent.

## 5. Branch B — `laderman23-orbit-min` (registered space and closure
argument)

**Question (fixed):** the exact minimum of
min over the 34,944 admitted data triples of the session-16 census of the
true three-stage total C(triple) = C(U) + C(V) + C(W_out), within the
frozen model and the +14 transposition accounting.

**Registered compute (complete census of exactly the undecided rows; no
sampling):**
1. Rows: from the hash-pinned session-16 artifact — per sigma class the
   mono row (125,125,125) at LB 59 and the 12 rows at LB 61; 39 rows
   total. For every row, its per-side d and floor-impossible flags are
   re-derived through this campaign's own d-counter/DFS and MUST equal
   the frozen values (assert; mismatch aborts).
2. Aux-1 push: every side of every one of the 39 rows has an impossible
   floor (frozen), so the complete next-level question per side is the
   aux-1 existential at T = d+1 over the pinned signed-pair-sum universe
   AU (closure lemma: at T = d+1 at most one gate is non-needed, so the
   auxiliary lies in AU; universes enumerated by two independent loop
   shapes that must agree). Dual instrument per instance: complete
   memoized subset-DFS AND kissat on `build_ext_cnf(tau + aux, T)`;
   100% agreement required; every infeasible instance carries a
   DRAT->LRAT proof accepted by BOTH pinned checkers; ledger per
   instance (CNF/LRAT sha256, rcs, seconds); checkpoint per instance
   (exact resume); ALL .lrat retained in the run dir.
   Instance count registered up front: 39 rows x 3 sides; |AU| = 408 at
   T = 15 (frozen pin for d = 14 sides) and |AU(T=16)| measured in-run
   for d = 15 sides (~450 by scaling, scouting estimate only).
3. Upper witnesses at the mono rows of sigma^1 and sigma^2 (LB 59 ->
   exact value): the frozen session-15 sigma^0 witnesses (16/16/16)
   transported by input-wire relabeling (the proved free-relabeling
   symmetry C(T(X)) = C(X), re-checked numerically in-run) to the
   sigma^1/2 mono sides, each accepted by exact expansion against that
   orientation's own 23 target rows; the output stage built by
   reverse-mode transposition of the transported Wfac-side circuit;
   assembled scheme verified END-TO-END by expansion over all 81
   monomials (9/9 outputs = the 3x3 matrix product) with the total
   recounted from the explicit gate lists. Expected 62 = 16+16+30; any
   other recount is a defect and aborts pre-verdict.
4. Verdict arithmetic: for each row, side value = d if floor achievable
   (none here), d+1 if the aux-1 census admits an auxiliary (witnessed),
   d+2 if it admits none; row total = sum of side values + 14; orbit-min
   = min over rows, cross-checked against the sigma^0 mono exact 62
   (frozen session 15) and the section-5.3 witnesses.

**Verdict rule:** if all censuses complete with 100% dual agreement and
100% certificate coverage: the exact orbit minimum is recorded
(FROZEN-CERTIFIED), expected 62; if any admission appears that yields a
row total < 62, the minimum is that row's exact value provided every
side entering it is exact (frozen floor-achievable, witnessed aux-1, or
pushed); if any side remains at an unresolved odd level (neither pushed
nor witnessed at some level < the incumbent), verdict FROZEN-INCONCLUSIVE
naming the rows. Budget: 16 CPU-h compute cap, 24 h wall cap; a bound
cap => FROZEN-INCONCLUSIVE with the exact decided prefix named (no
partial negatives, ever).

## 6. What this campaign may NEVER claim

No universal no-54 (or no-55) claim about rank-23 3x3 multiplication;
nothing about decompositions outside the named candidate list and the
five frozen triples; nothing about non-ternary alphabets, GL(3,Q)/GL(3,Z)
beyond the registered families, anti-cyclic BA swaps, actions outside
the proven automorphism family, or multi-auxiliary spaces (a Branch-B
negative at aux-1 leaves aux-2+ untouched and says so); Branch A's
negative is a statement about D6's registered orientation space only.
A Branch-B verdict of 62 is a statement about the `laderman23` orbit
inside the frozen model; it is NOT a statement that 55 is optimal, and
Laderman's 62 stays above the record 55.

## 7. Budget summary (fixed now)

| phase | cap |
|---|---|
| A: anchors + controls | 0.5 CPU-h |
| B: source audit (fetches, transcription, Brent x2, equivalence battery) | 2.0 h wall |
| Branch B census (39 rows x 3 sides x |AU|, dual instrument + certificates) | 16 CPU-h |
| Branch B witnesses + end-to-end | 1 CPU-h |
| postcompute audit | 1 h wall |
| Branch A census (if it fires) | 14 CPU-h |
| total wall | 24 h |

Measured rates (session 15: 1,224 dual instances + certificates in
~7.5 CPU-min; session 16: 8,208 dual pair decisions + 7,029 certificates
in ~1.2 h wall) make Branch B's ~51,000 instances the long pole at a
projected 5-8 CPU-h. If any cap binds: FROZEN-INCONCLUSIVE with the
exact decided prefix named and the remainder explicitly out of scope.

## 8. Reproduction (fixed)

Run dir is self-contained after freeze (scripts as run, prereg copy,
tools, ledgers). Census rerun:
`<venv>/bin/python orbit_min_run.py --run .`; audit rerun:
`<venv>/bin/python orbit_min_run.py --run . --phase audit`; postcompute:
`<venv>/bin/python postcompute_audit.py`.

## 9. Terminal verdicts (exactly one)

- **FROZEN-CERTIFIED** (Branch B, expected): orbit-min = 62 exact (or a
  lower exact value fully witnessed), all controls pass, 100% dual
  agreement + certificate coverage on the complete registered census.
- **FROZEN-NEGATIVE** (Branch A only): the sixth-source registered
  orientation space contains no <= 54-certifiable triple; histogram and
  minimum as reported.
- **FROZEN-INCONCLUSIVE**: any cap binds, any instrument disagreement or
  certificate gap survives retry, any control fails after an amendment
  cycle — with the exact reason and decided prefix named.
- LIVE branch: any certified LB <= 54 (Branch A) or any witnessed total
  < 59 (Branch B, impossible short of an instrument defect) halts before
  recording and escalates to Main via hub.
- Control failure: abort with NO verdict; amendment only with both
  prereg hashes recorded before any amended rerun.
