# Theoretical CS / information-theory open-problem attack repository

**Start here:** [`RESULTS.md`](RESULTS.md) is the authoritative current-results
index. This README is the repository map. Provisional claims and candidate
results live in a non-authoritative `<target>/scratch/` subfolder until verified.

## Mission

Attack open problems in **theoretical computer science and information
theory** — complexity, algorithms, coding, proof complexity, cryptographic
hardness — with the same verification discipline as the sibling `../math/` and
`../physics/` repositories.

The selection bias is explicit and load-bearing: prefer problems whose progress
is expressible as a **machine-checkable finite certificate** (DRAT/LRAT, VeriPB,
exact rational LP/SDP dual, exhaustive census with a completeness argument) on a
single workstation, over problems whose only currency is asymptotic argument or
floating-point benchmark. A solo repo cannot win a scaling race; it can own the
certification niche that funded groups have no incentive to build.

## Tracking contract (SSOT)

One fact, one owner. Four layers, no duplication:

| layer | owner | contents |
|---|---|---|
| authoritative results index | [`RESULTS.md`](RESULTS.md) | one row per target: status, headline result, evidence pointer; campaign ladder; adopted next targets |
| chronological ledger | [`PROGRESS.md`](PROGRESS.md) | dated, newest-first session entries; every claim carries its verification; retractions inline |
| current state per target | `<target>/README.md` | what is PROVED / CONDITIONAL / NUMERICAL right now, gates, file inventory, how to run |
| run artifacts | `<target>/campaigns/<UTC-timestamp>_<uuid>_<code-hash12>/` | immutable frozen snapshots + committed results; inventoried in `<target>/campaigns/README.md` |
| generated compact state, per target | `<target>/state.json` | owner-maintained semantic fields (`problem`, `current_gate`, `headline`, `next_action`); mechanics refreshed by `scripts/campaign.py state refresh` |

Generated, not narrative: when the one-line summary in `<target>/state.json` and an
owner document above disagree, the owner document wins — fix [`RESULTS.md`](RESULTS.md),
`<target>/README.md` (its append-only rule is unaffected), or [`PROGRESS.md`](PROGRESS.md)
first, then re-run `scripts/campaign.py state refresh` so the generated copy follows.

Naming rules:

* **One problem, one top-level folder**, short lowercase name. No dated
  top-level folders, no second folder for the same problem.
* Campaign directories are producer-generated (`<UTC>Z_<uuid>_<hash>`); never
  hand-created, never edited after launch. Disposable rehearsals go under
  `<target>/campaigns-smoke/`.
* Update the ledger and the owning README **in the same session** as the work;
  a result that is not in `PROGRESS.md` does not exist.
* Unverified/candidate material lives under `<target>/scratch/`
  (non-authoritative); promotion to a README requires the evidence labels below.
* `docs/scan-raw/` holds raw subagent scan output verbatim. It is
  **non-authoritative** and is never cited as evidence; the distilled,
  owner-checked scan document is.
* **Target READMEs are append-only below the contract.** A target's
  `README.md` serves two roles at once: it is the authoritative gate contract
  (claim, gate ladder, pre-campaign requirements, disjointness) and it is the
  per-target current-state document in the table above. An owning agent
  **appends** a `## Current state (agent <Name>, <date>)` section carrying
  verdicts, frozen campaign paths, evidence labels, honest gaps and how-to-run.
  It never deletes, replaces or truncates the contract. If a gate's premise is
  wrong, a correction note goes *beside* it, never in place of it.
  Established 2026-08-30 after agent `Mm3` replaced `mm3/README.md` with a
  campaign summary, deleting the gate ladder; the ladder was restored and the
  agent's content retained verbatim beneath it.
* **Corrections carry attribution, in both directions.** When an owning agent
  shows that something in a gate ladder is wrong, the owner fixes it inline and
  names the agent. Four such corrections landed on 2026-08-30: `omega/`'s
  "published constants are all float-only" line (wrong — the three predecessor
  rungs released parameters on OSF), `rs-pe3d/`'s "$12\mid q-1$" instance
  condition (wrong — the requirement is three pairwise coprime divisors of
  $q-1$), `mceliece/`'s instance ladder (infeasible — violated $k=n-mt>0$), and
  `kg/`'s gate-A runtime estimate (optimistic). Conversely, when an agent
  misattributes its own premise to the owner, that is corrected too:
  `rs-pe3d/`'s false "Anchor 1" equality was its own derived unit premise, not
  a README claim.

## Targets

The initial seven targets opened 2026-08-29 from the ranked shortlist in
[`docs/CS_FRONTIER_SCAN_2026-08-29.md`](docs/CS_FRONTIER_SCAN_2026-08-29.md).
Per-target state, gates and pre-campaign requirements live in each
`README.md`; verdicts are collected in [`RESULTS.md`](RESULTS.md).

**The third column below is the problem framing when each target was opened. It is
deliberately not maintained.** Current status is owned by [`RESULTS.md`](RESULTS.md) and by each
`<target>/README.md`; duplicating it here would create a second authority for one fact, which the
tracking contract above forbids. So `BENCHMARK` in this table means "this is what the target was
opened to do", never "this is where the target is now" — several have since passed gates A, B
and C. Read the column for the *problem*, follow the link for the *state*.

| directory | problem | framing at opening — NOT current status |
|---|---|---|
| [`mceliece/`](mceliece/README.md) | **Classic McEliece hold-out "waterfall" dispute** — is the derivative-flag collapse at $c=2t+3$ real on binary Goppa? | BENCHMARK. Live four-paper ePrint dispute (2026/1747 → 2026/1810 → 2026/1786 → v2), all abstracts owner-read. Gate A is a nine-instance exact $\mathbb F_{2^m}$ nullity census, under an hour, decisive either way. |
| [`kg/`](kg/README.md) | **Grothendieck constant $K_G$** | BENCHMARK. Triple-arrival target. $6\pi/11\le K_G\le\pi/(2\log(1+\sqrt2))-3.47\times10^{-4}$ (arXiv:2608.11158, owner-read); the tenths-digit corollary is already re-derived here in Arb. Both bounds reduce to finite interval-arithmetic inequalities; gate C attempts a certified improvement. |
| [`omega/`](omega/README.md) | **matrix multiplication exponent $\omega$** | BENCHMARK. $\omega<2.371177$ (arXiv:2608.16884, owner-read) comes from an ML+AlphaEvolve numerical optimum with no published rigorous enclosure. Gate A asks whether the record is rigorously established as stated. |
| [`delcap/`](delcap/README.md) | **binary deletion channel capacity bounds** | BENCHMARK. Whole published constant chain (0.1185, 0.4143, 0.3745, 0.1221, 0.3578, small-$d$ expansions) is floating-point Blahut–Arimoto; none ever enclosed rigorously. Rubinstein–Con owner-read. |
| [`oct-rank/`](oct-rank/README.md) | **real tensor rank of octonion multiplication** | BENCHMARK. $18\le\mathrm R_{\mathbb R}(T_{\mathbb O})\le25$ (arXiv:2608.16649, owner-read), narrowed from 15–30 twelve days ago. Gate B attacks rank 24 by Krawczyk certification; lower-bound search is explicitly ruled unsound. |
| [`rs-pe3d/`](rs-pe3d/README.md) | **high-dimensional product expansion for Reed–Solomon tensor codes** | BENCHMARK. ECCC TR26-150 Conjecture 4.2 (owner-read abstract); lowest duplication risk in the scan — the only existing attempt is an AI proof the authors state is unverified. Near-linear private PCPs depend on it. |
| [`mm3/`](mm3/README.md) | **additive complexity of rank-23 $3\times3$ matmul** | BENCHMARK. 55 additions (arXiv:2607.28676, owner-read); the "provably optimal for this fixed orientation" claim is ILP/SAT-decidable, and 54 would break the record. |
| [`xor/`](xor/README.md) | **Exact XOR synthesis for AES MixColumns** (opened 2026-09-06) | Assay of local headroom in the public 88-XOR circuit; exact GF(2) functional replay and independently bounded local resynthesis, not a private benchmark or an AI-superiority claim. |

## Target-selection surveys

* [`docs/CS_FRONTIER_SCAN_2026-08-29.md`](docs/CS_FRONTIER_SCAN_2026-08-29.md)
  — six-scout parallel frontier scan (meta-complexity, fine-grained complexity,
  information/coding theory, proof complexity + certified search, cryptographic
  hardness, learning/communication complexity), with the ranked candidate
  shortlist and the first falsifiable gate per candidate.
* `docs/scan-raw/` — the six raw scout reports, verbatim, non-authoritative.

## Discipline

Rules 1–7 inherited from `../math/README.md` without modification; 8–13 are TCS-specific;
17b. **A count is not evidence until its pattern is anchored.**
    Added 2026-08-30. Bare substring counts inflate reference counts with ordinary words —
    `ns` matched 140 lines of *constraints*/*functions*, `dist` matched *distance*/*distinct*;
    both were 0 when anchored to `` `t/` ``. Anchor the pattern, then report the count. And
    before calling a difference a defect, check whether the difference is **declared**: an
    audit that flags a documented, deliberate divergence is measuring the wrong thing. Both
    failures occurred in one owner check on the same day, and the earlier observation that
    opened that thread had used the same unanchored method — right by luck, wrong in detail.

17c. **Verify provenance before attribution. A correct measurement of the wrong object is
    still an attribution error — and a clean number is what disguises it.**
    Added 2026-08-30. The owner audited a vendored external Lean checkout competently (492
    files, zero `axiom` declarations, zero `sorry` in the load-bearing chain, `sorry`s confined
    to unimported challenge files) and then added it to a sibling repository's Targets table —
    which would have attributed an outside unconditional formalization of $H_1\le246$ to that
    repository. The directory was **already documented** as third-party, in a section the owner
    had written hours earlier, saying explicitly *"do not cite as a result of this repository."*
    Before attributing an artifact: check for a vendored `.git/`, an upstream URL, a
    `lean-toolchain` or lockfile you did not author, and read the existing documentation of the
    thing you are about to describe. Rules 14-17b guard the instrument and the domain; this one
    guards the **object**, and it is the one that turns good arithmetic into a false claim.

**14–17c were added 2026-08-30 and are the ones this repository learned the hard way** — each
records a real failure from that session, two of them in the owner's own reasoning.

1. Read primaries first-hand; quote verbatim. Secondary summaries are not evidence.
2. Every quantitative claim is machine-checked by a script in this repo.
3. Distinguish PROVED / CONDITIONAL / NUMERICAL at every step, always.
4. Verify by hand anything an optimiser reports before believing it.
5. Record retractions inline rather than silently editing claims away.
6. A theorem is not proved because its asserts pass — asserts check arithmetic at
   the points you chose.
7. **Never promote a sampled parameter check to a universal claim.** Where the
   quantity is monotone or affine in the parameter, solve for the boundary and
   assert the equivalence. Sampling a range is evidence, never a quantifier.

TCS-specific additions:

8. **Every complexity claim states its model exactly**: machine model,
   uniformity, alphabet, worst-case vs average-case, and — for an asymptotic
   claim — the regime and whether constants are hidden.
9. **Finite-n verification is never an asymptotic claim.** "Checked for all
   $n\le N$" is a finite theorem; extrapolating it is a conjecture, and a fitted
   exponent is not a bound. State $N$, state the conjecture separately.
10. **Every implemented algorithm's claimed complexity is measured, not
    assumed**: report an operation counter, not wall-clock alone; any wall-clock
    number names CPU, thread count, and toolchain version.
11. **Search results carry their certificate.** UNSAT ⇒ a DRAT/LRAT or VeriPB
    proof accepted by a pinned checker, with the checker version recorded;
    exhaustive census ⇒ an explicit canonical-form/orbit completeness argument;
    heuristic search ⇒ `COMPUTATIONAL-EVIDENCE` only, never a bound.
12. **Conditional results keep their hypothesis in every restatement.** Name the
    conjecture (SETH, ETH, 3SUM, LWE, ...); a result that drops its assumption
    in the abstract is a different, false result.
13. **Reproduction is not re-derivation.** Tag a reproduced published number
    `[REPRODUCED]` (same inputs, same stated pipeline) versus `[DERIVED]`
    (independent implementation from the definition). They are different claims.
14. **THE INSTRUMENT MUST BE ABLE TO ESTABLISH THE CLAIM, NOT MERELY BE CONSISTENT
    WITH IT.** Added 2026-08-30 after this failure shape appeared **four times in one
    session**, twice in the owner's own reasoning. Before believing any check, ask what
    it *could not* have detected:
    * **Shared component.** A validation that shares machinery with the thing it
      validates cannot see a defect in that machinery. `kg/`'s finite-difference audit
      checked a wrong function against the analytic derivative of *that same wrong
      function* and passed by construction. Corollary: $k$ solvers over **one** encoding
      are $k$ checks of the encoding's *consequences*, not of the encoding.
    * **Partial coverage stated as complete.** Finitely many rejected directions say
      nothing about a positive-dimensional space; an *inscribed* region is not its
      circumscribing one. A no-go over a set of directions is worth exactly the **rank
      argument** behind it — if the rank equals the coordinate count the feasible set is a
      point and no probe is needed; if it is less, no finite number of probes settles
      anything. There is no middle case, and *"the mechanism looks structural"* is not a
      substitute for the count.
    * **Representation artifacts.** A round number in a scaling factor is a basis artifact
      until proven geometry. An inscribed cube in an unnormalized basis certifies a domain
      shrunk by roughly the dimension.
    * **Tolerance versus quantity.** A solver whose feasibility tolerance equals the
      magnitude being constrained certifies nothing: HiGHS at its default $10^{-7}$, asked
      to respect bounds of size $10^{-7}$, returned a point violating them by 100% with a
      plausible objective 14% off. **Non-dimensionalize**; never merely tighten tolerances.
      And **re-verify every solver's returned point against its own constraints in
      independent arithmetic** — a solver's reported objective and its reported feasibility
      are both claims, not evidence.
15. **Name the certified domain in every certified sentence.** "Certified" without a named
    domain is how overclaims enter. If the domain is an inner approximation, state the
    uncovered remainder explicitly. Added 2026-08-30 after a local-optimality claim was
    written and retracted the same day.
16. **Domain-shopping is prohibited.** A domain, box, search space, or candidate set is
    fixed *before* the run and reported with every outcome including misses. Moving or
    shrinking it after seeing the result — until the bound closes, or until the target
    value is hit — is the numerology failure mode at the level of search spaces, and it is
    prohibited whether or not a citation is attached. If a region must change, the new
    region is a **separate** pre-registered statement reported alongside the original and
    its outcome.
17a. **Sources without a dated API are tracked by CONTENT HASH, not by a listing date.**
    Added 2026-08-30. Rule 14's companion for provenance: the arXiv API gives a per-version
    `updated` field, so arXiv sources can be re-checked by date. IACR ePrint and ECCC do not
    expose an equivalent, and a listing date can lag or be absent. So **cache the artifact and
    compare bytes on every re-check** — size, hash, and text line count. This is strictly
    stronger than a date, because it detects a silent in-place revision that no date would
    reveal. It caught a real one the day it was adopted: ePrint 2026/1747 was revised in place
    on 2026-08-30, detected as 245,667 bytes / md5 `cd044d4c…` becoming 264,802 bytes / md5
    `f27d10be…`, text growing 1633 to 1835 lines. **Keep the superseded copy alongside the new
    one under a VERSION-STAMPED filename** (e.g. `pe1747_v3_2026-08-30.pdf` beside `pe1747.pdf`),
    never overwriting a cached primary — so that the **diff stays reproducible later**, not merely
    the detection. Refinement contributed by agent `Mceliece` from applying the rule.)
    **2026-08-31 clarification after the rule itself was misapplied:** the unqualified
    canonical path (for example `pe1747.pdf`) MUST contain the current fetched bytes.
    Every displaced content version gets its own source-version or timestamp-stamped path;
    the stamp is archival identity, not a synonym for “current.” `mceliece/` had these roles
    reversed for `1747`, and the error survived because both bytes were retained but the
    README named the wrong one canonical. “Never overwriting” above is a **byte-retention**
    invariant, not a frozen-path invariant: first archive and hash the displaced bytes, then
    replace the canonical path with the current bytes. Both versions remain addressable and
    checksum-verifiable. At campaign startup, also inspect any source-provided
    version history: this caught four later `1786` revisions after a byte-stable cached
    `20260828:000118` copy. Record SHA-256, byte count, and source revision timestamp where
    available. Retain every materially used version; never repair provenance by deletion.
17. **Owner-supplied figures are inputs, not evidence.** A number handed down by the
    session owner carries no more authority than one from a paper: re-derive it before any
    claim rests on it. Added 2026-08-30 after two owner figures propagated through agents
    unchecked in a single session — one a structural hypothesis, one a rank count that was
    wrong.
17d. **A boundary is per-instrument, and an absence probe must name what it bounds.**
    Two distinct failures, both found 2026-08-31 in this repository's own sweeps, both
    invisible to arithmetic checking:
    * **One boundary per category, printed.** The 2026-08-30 sweep printed a single
      `max published` (`2026-08-27T17:53:45Z`) for a **two-category** sweep. That value is
      `cs.CC`'s stratum maximum — re-verified exactly — while `cs.IT`'s 08-27 maximum is
      `2026-08-27T20:14:37Z`. Three `cs.IT` submissions lived in the 2 h 21 min shadow
      (`2608.27565`, `2608.27635`, `2608.27682`), and the following window's query then
      excluded the whole 08-27 day by its own lower bound, so **two instruments missed the
      same items for two different reasons**. A sweep over $k$ sources or categories owes
      $k$ printed boundaries; a maximum generalized across categories is a coverage defect,
      not a summary. This is rule 17b's anchoring requirement applied to the boundary itself.
    * **Soft-404: status is not content.** `eccc.weizmann.ac.il/report/2026/163` returns
      **HTTP 200** carrying a `404 PAGE NOT FOUND` body. A status-based absence probe
      concludes the report EXISTS; only a content discriminator settles it (TR26-162's page
      carries `TR26-\d+` identifiers, TR26-163's carries none). State which layer the probe
      reads — transport status or document content — in the absence claim itself.
17e. **A rate is a claim about an accounting source; name it.** Added 2026-09-01
    after a 40x phantom CPU clamp cost three agents and the owner real hours. On this
    workstation `ps -o time` / `%cpu` for framework-Python (`Python.app`) children
    quantizes and defers CPU-time accounting: repeated 25 s windows read `0->1 s`
    while `/usr/bin/sample` showed the same thread at `2287/2287 ms` inside a numpy
    ufunc, i.e. a full core, with throughput independently confirming it. Four
    unrelated processes all reporting exactly `2.4%` was the tell that the
    *instrument* was constant, not the load. Trustworthy sources, in order:
    in-process `resource.getrusage(...).ru_utime` or `time.process_time()`; then
    `/usr/bin/sample <pid> N` thread time and state; then throughput against a
    calibration run. The owner's own prescribed discriminator (`delta_cpu/delta_wall`
    from `ps`) was the broken one, and the owner's A/B tests were right only because
    they happened to use in-process accounting. Corollary, and the reason this is a
    rule: **never restructure a campaign around an unexplained performance
    observation** - get a stack first, and retain it before killing the process, or
    the evidence dies with it (it did, twice, that day).

### Verify before record — the mechanism, not the intention

Rules 2, 17 and 17c are enforced at the moment a number enters a ledger, which is
where they historically failed: under batch load the owner records a returned
figure instead of re-deriving it. The standing protocol, unchanged since the
seven-target launch, is that **this session (Main) is the single writer** of
`README.md`, `RESULTS.md`, `PROGRESS.md` and `docs/`; an agent writes only inside
its own `cs/<target>/`, and `docs/scan-raw/` output is non-authoritative because a
scout tag is subagent-verified, not owner-verified.

Recording is gated by [`tools/verify_campaign.py`](tools/verify_campaign.py), an
owner-only reader of finished artifacts (no campaign runner may import it). Given
a campaign directory it re-derives, from bytes on disk: every checksum-ledger
entry; every newline-inclusive row hash against the per-line ledger, with the
row/ledger key bijection and ordinal order; the declared source hashes against the
files that carry them; the outward decimal ↔ archived binary-rational identity and
the enclosure of each archived Arb ball by the outward endpoints and width; each
presentation CSV as an order-independent row-wise projection of a canonical
`.jsonl`, reporting **every** row field a column could have come from; and the
sums, maxima and histograms of every leaf shared by all rows. `--claims FILE`
diffs stated numbers against those derived values and exits non-zero on the first
disagreement, so a ledger sentence is written only after the numbers in it have
been reproduced from the frozen files.

Two negative controls were run when it was introduced, both on copies: perturbed
claims (an off-by-one count, a wrong row total, a wrong file hash) were each
reported with the artifact value, and a single flipped byte inside a canonical row
was caught independently by the file ledger, the per-line hash, and the CSV
projection. Unlisted campaign files are always reported, so an artifact that no
ledger covers cannot quietly back a claim. A campaign whose runner never promised
sorted-compact canonical rows is reported, not failed — rule 17b: a declared
difference is not a defect.

## Evidence labels

Literal, as in `../math/RESULTS.md`:

* **MACHINE-VERIFIED** — a named executable fact: a script in this repo, or a
  pinned external checker, accepts it.
* **HUMAN-AUDITED** — ordinary checked mathematics, line by line.
* **COMPUTATIONAL-EVIDENCE** — non-proof numerical evidence.
* **CITED-DEPENDENCY** — imported from a primary, with the primary read.
* **OPEN** — unresolved.
* **FAILED** — the wording or inference is not established (retraction).

Plus the provenance tags `[REPRODUCED]`, `[DERIVED]`, `[REPORTED]` (not read
first-hand), `[INFERENCE]` (inferred, not read), `UNVERIFIABLE` (citation could
not be fetched).

## Resource policy

- One low-priority process (`nice -n 10`) at a time.
- Pin OMP, OpenBLAS, MKL, vecLib, and NumExpr to one thread.
- Bounded searches with explicit stop conditions; no daemonized CPU loops.
- Randomized search: fixed seed per run, seed recorded in the campaign
  inventory; repeated-seed policy defined before the first run.
- Every SAT/PB run records solver version, encoding hash, and checker version.
- Run prior-art review before announcing novelty. Keep failed routes and
  corrections in the chronological ledger.

## Layout

| Path | Contents |
|---|---|
| `RESULTS.md` | Authoritative per-target results index. |
| `PROGRESS.md` | Chronological session ledger (newest first). |
| `<target>/` | One folder per problem, when open. |
| `<target>/README.md` | Per-target authoritative state and gates. |
| `<target>/campaigns/` | Frozen run snapshots, producer-generated names. |
| `<target>/scratch/` | Non-authoritative candidate material. |
| `docs/` | Frontier scans and background reading. |
| `docs/scan-raw/` | Raw subagent scan output, non-authoritative. |

## Evidence discipline

1. Every claim carries a source: arXiv ID, ECCC report number, DOI, or an
   in-repo derivation with committed campaign artifacts.
2. Numbers copied from papers name the figure/table and the identifier.
3. Not personally verified against a primary source → tag `[REPORTED]`.
4. Inferred rather than read → tag `[INFERENCE]`.
5. Conflicting values are recorded as conflicts, not silently resolved.
6. Numerical candidates remain `NUMERICAL`; only exact/verifiable evidence
   (proof certificate accepted by a pinned checker, exact rational arithmetic,
   Lean proof, machine-checked replay of a full campaign) promotes to `PROVED`.
7. **A no-go theorem is worth exactly its hypotheses.** Transcribe them; the
   surviving hypotheses are where a construction can live.

## Local-machine resource contract (owner instruction, 2026-09-04)

The owner requires total CPU below 50% and preservation of local disk
headroom. Launch numerical work through [`tools/resource_guard.py`](tools/resource_guard.py)
under the process supervisor, never through an unbounded background shell.
Only audited single-threaded commands that wait for their children may use
this wrapper. Its constants own the operational thresholds: two concurrent
slots, CPU pause/resume margin below the owner's ceiling, a 100 GiB disk
reserve, bounded file/output growth, and a sampled process-group RSS stop.
Every job also has a finite wall ceiling. Resource stops preserve partial
evidence and are never successful scientific outcomes.

The guard signals only the process group it created. It cannot impose a
global quota on unrelated applications or prevent their instantaneous CPU
spikes; samples and pause/stop events are the verification evidence. Do not
bypass it with scheduler boosts, additional workers, large artifact copies,
or unregistered retries. Keep scope-specific resource amendments beside each
live campaign, preserve failed attempts, and do not delete existing files to
recover space. [`tools/test_resource_guard.py`](tools/test_resource_guard.py)
exercises the stop/admission boundaries using small real processes without
generating CPU or memory pressure.

## Autonomous research supervisor (owner-authorized, 2026-09-05)

A durable macOS user LaunchAgent (`local.jinleic.cs-research`) dispatches one
bounded research worker at a time against the queue in
`~/Library/Application Support/cs-research`. Code:
[`tools/research_supervisor.py`](tools/research_supervisor.py) (service,
recovery, promotion), [`tools/research_queue.py`](tools/research_queue.py)
(durable SQLite control/tasks/compute leases) and
[`tools/research_tools.py`](tools/research_tools.py) (the only tools a worker
has: bounded read/list/search, draft writes, preregistration, guarded compute,
one structured finish).

Control, from the workspace root:

```
python3 -B cs/tools/research_supervisor.py status
python3 -B cs/tools/research_supervisor.py stop   --reason "..."
python3 -B cs/tools/research_supervisor.py resume --reason "..."
```

`stop` is sticky: it disables the store, bumps the control generation, stops
owned worker and compute groups by kernel birth identity, and survives
restarts, reinstallation and crashes. Only an explicit `resume` re-enables
dispatch. Nothing is deleted on stop.

Invariants the machinery enforces, not the model's cooperation:

- **Preregister before compute.** A statement plus every source is committed
  locally, copied into the producer run and hashed; sealed bytes are immutable
  and unreadable-if-modified thereafter. An interruption after producer init is
  recovered from a durable intent journal, never re-minted.
- **Guarded, sandboxed compute only.** Sealed scripts run under
  `tools/resource_guard.py` (hash-pinned) plus a deny-default `sandbox-exec`
  profile: writes confined to the run dir with sealed sources, producer
  metadata and importable code denied; no network, no fork, no hardlinks. Each
  compute is probed before launch and refuses if any boundary is open.
- **No unchanged reruns.** Identical sealed source plus logical argv (inputs by
  hash) can never execute twice; only a never-started resource wait is reusable.
- **Independent verification before promotion.** `FROZEN-CERTIFIED` requires a
  candidate/complete primary with a successful receipt, plus a separate review
  job with its own preregistration, its own successful receipt, a `verified`
  outcome and a binding to the primary commit. An incomplete review retains
  the primary for a later review, rather than inventing a scientific verdict.
- **Honest recovery.** No structured finish means no infrastructure-authored
  scientific outcome. Never-started attempts defer with their sealed registration;
  uncertain or executed-but-unreported primaries go to independent review and
  cannot be certified without a primary finish. Receipt execution status stays
  authoritative; storage-limit observations are separate facts.
- **Bounded execution and pauses.** A compute has at most 600 s charged execution,
  a prelaunch admission wait no longer than its requested wall, and at most
  900 s cumulative pause time. The watchdog also bounds total real elapsed time.
  A worker has 990 s charged time (600 compute + 90 guard grace + 300 turn reserve)
  and at most 1500 s cumulative pauses. Pauses do not consume execution time;
  finite pause ceilings retain partial evidence and release owned processes.
  Live remaining compute budgets include both limits and are rechecked after probes.
  Other ceilings remain: 64 MiB growth per compute, 96 MiB cumulative evidence,
  2 GiB sampled group RSS, 100 GiB disk reserve, and capped logs. Resource admission
  waits retry after 30 s without increasing failure backoff; infrastructure/failure
  backoff is exponential, capped at 900 s.
- **Owner-state preservation.** Closeouts update only `autonomy_latest` and
  producer-derived fields, carrying all existing owner fields through state import.
  Human-owned inventories are preserved and refresh skips are reported explicitly.

Verification: `tools/test_research_supervisor.py` and
`tools/test_resource_guard.py` exercise real process suspension/resumption,
finite pause expiry, guarded sandboxed compute, receipt/queue consistency,
restart recovery, retained registrations, independent-review dispatch and sticky stop.
Host-load boundaries use scripted samples; the tests do not generate host load.
Deployment evidence is retained under
`.autonomy/supervisor-repair-20260905T074741Z/`.
