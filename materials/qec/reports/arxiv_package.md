# arXiv submission kit — paper_pbb_nogo

**DO NOT POST YET.** This file is the submission staging document only; no arXiv
submission has been made. Compiled: 2026-08-24 against `paper_pbb_nogo.md`
(current revision, abstract with **two** boxed results; verified verbatim, see
§9) and `paper_pbb_nogo.pdf` (pdflatex: 18 pages, 716,216 bytes).

---

## 1. Title and abstract

**Title field (verbatim):**

    Perturbing a bivariate-bicycle code cannot buy distance for free

**Subtitle** (second line of the paper masthead; not part of the arXiv title field —
fold into the abstract's first sentence only if a combined one-liner is wanted):

> An exact, perturbation-independent rate–distance obstruction for perturbed bivariate-bicycle (PBB) codes, with machine-checked certificates

**Abstract field (verbatim from `reports/paper_pbb_nogo.md`, lines 13-63,
`## Abstract` body, hard-wrapped exactly as in the source; contains the TWO boxed
results. Programmatic verbatim check: PASS — see §9):**

```text
Perturbed bivariate-bicycle (PBB) codes are non-CSS stabilizer codes obtained from a CSS
bivariate-bicycle (BB) parent $P$ by adding a $Z$-type perturbation $[C\;D]$ to the
$X$-check block. They have been proposed as a route past the CSS BB rate–distance
envelope, on the strength of a published catalogue of $368$ codes including a non-CSS
$[[144,12,12]]$ and a $[[360,12,\le 24]]$ candidate.

We prove that the perturbation degree of freedom is not free. Let
$\Delta=\{\lambda[C\;D] : \lambda[A;B]=0\}$ be the *dressing space* of the perturbation.
We show (i) $\Delta$ is an $R$-submodule of the $Z$-sector for
$R=\mathbb{F}_2[x,y]/(x^{\ell}-1,y^{m}-1)$, (ii) the PBB dimension obeys the exact
identity $k_Q = k_P - \dim\bar\Delta$, and (iii) any minimum-weight $Z$-logical of the
parent that is *not* absorbed by $\Delta$ survives as a logical of $Q$, certifying
$d_Q \le d_Z(P)$.

Because $\Delta$ is a submodule, absorbing one minimum-weight parent logical absorbs its
entire translation orbit. This upgrades the per-perturbation criterion to a
**parent-level no-go theorem**: with
$T(P)=\dim\big(M(P)+S_Z\big)/S_Z$, where $M(P)$ is the span of the translation orbits of
*all* minimum-weight $Z$-logicals of $P$,

$$\boxed{\;d_Q > d_Z(P)\ \Longrightarrow\ k_Q \le k_P - T(P)\quad\text{for \emph{every} perturbation } [C\;D]\;}$$

$T(P)$ depends only on the parent and is computed *exactly* — not bounded — by a SAT
enumeration that terminates in a certified UNSAT. When $T(P)=k_P$ the entire family is
closed: no PBB over $P$ retains a single logical qubit while exceeding $d_Z(P)$.

We compute $T$ exactly for $134$ of the $202$ distinct parents in the published
catalogue — exhaustively for every parent at $n\le 144$. **$61$ parents are family-closed**, including the Gross code
$A=x^3+y+y^2,\;B=y^3+x+x^2$ itself ($T=12=k_P$) and all $11$ distinct catalogue
$[[144,12,12]]$ parents. This settles $249$ of the $368$ catalogue rows with no per-row
search at all, and it forecloses the headline direction of the construction over the
most important parent in the family.

The bound is not merely valid but *exact*.  A second elementary identity — the
dressing space is the image of the left kernel of $[A\;B]$, whose dimension is
$k_P/2$ on every BB parent — caps $\dim\bar\Delta\le k_P/2$.  Whenever
$T(P)\ge k_P/2$, which holds on $133$ of the $134$ certified parents, the two
bounds **sandwich**:

$$\boxed{\;d_Q > d_Z(P)\ \text{and}\ T(P)\ge k_P/2\quad\Longrightarrow\quad
k_Q = k_P - T(P)\ \text{exactly}\;}$$

Perturbation does not just pay at least $T$ logical qubits for distance; on
every parent in the theorem's class it pays *precisely* $T$, no more and no
less.  All $7$ independently certified distance-increasing perturbations sit on
the law exactly (slack $0$, hypotheses replay-certified), and a controlled
small-lattice probe — $64$ sampled parents across four lattices, $26{,}898$
valid $\delta>0$ perturbations of $Z$-support $\le 4$, distances decided by
exact symplectic meet-in-the-middle — found $496$ further
strict increases, **every one on the law, zero above it**: instances we now
know were theorem-forced rather than lucky.
```

*arXiv LaTeX note:* the abstract renders under arXiv's math support as-is
(`$...$` + `$$\boxed{...}$$` are fine; there is one display-math block per boxed
result, and `\text`/`\emph` usage is standard). No cut needed: 3,158 chars ≈ 500
words — inside arXiv's abstract length limit.

---

## 2. Categories

| Slot | Category | Rationale |
|---|---|---|
| **Primary** | `quant-ph` | Quantum error-correcting codes community: BB/Gross code (arXiv:2308.07915), the PBB catalogue (arXiv:2606.02418), and all cited QEC literature are quant-ph; the claims are physical-code (rate–distance, fault tolerance) claims, so quant-ph indexing reaches the intended readership and review pool. |
| Cross-list | `cs.IT` | The result is a rate–distance obstruction for a stabilizer-code family proved by exact GF(2) linear algebra and SAT enumeration; cs.IT is the home of information-theoretic coding bounds and reaches the algebraic-coding audience that consumes envelope/no-go results. |
| Cross-list | `math.CO` | The machinery is combinatorial: modules over the ring $F_2[x,y]/(x^\ell-1,y^m-1)$, translation orbits of minimum-weight logicals, and an exact counting theorem; math.CO covers combinatorial matrix theory / designs-and-codes work of exactly this character. |

No `math.QA`/`cs.ET` cross-list: the paper makes no operator-algebra or
hardware-architecture claim.

---

## 3. Comments line

Suggested `Comments:` field (measured against the compiled PDF, not the md):

> 18 pages, 7 tables, 0 figures; two boxed main results; machine-checked with
> exact SAT and rooted connected-cluster certificates plus GF(2) replay;
> full provenance ledger and per-claim reproducers in the ancillary bundle

Measured, not guessed: pdflatex line
`Output written on paper_pbb_nogo.pdf (18 pages, 716216 bytes)`; 6 table floats
plus 1 longtable, 6 tabular environments plus 1 longtable, and 0 figures.

---

## 4. Ancillary files

Upload the tex sources normally plus these as ancillary files (arXiv "ancillary files"
mechanism). Paths relative to repo root `math/qec/`; sizes as of 2026-08-24.

| File | Size | One-line description |
|---|---|---|
| `reports/paper_pbb_nogo.tex` | 71,077 B | LaTeX source compiled from the md (arXiv requires source; ship as the primary tex). |
| `reports/paper_pbb_nogo.bib` | 3,173 B | Bibliography (compiled `paper_pbb_nogo.bbl` produced from it; arXiv will run its own bibtex pass). |
| `reports/paper_pbb_nogo.pdf` | 716,216 B | Reference build of the paper (18 pages) for side-by-side preview; arXiv rebuilds from tex. |
| `results/certificates/exp056_wm_162_8_14_distance.json` | 127,917 B | EXP-056: exact $[[162,8,14]]$, 20 replayed class-orbit UNSATs plus weight-14 witness and duality. |
| `results/certificates/exp057_170_16_10_distance.json` | 4,516 B | EXP-057: replayed exact $[[170,16,10]]$. |
| `results/certificates/exp058_186_10_14_distance.json` | 5,572 B | EXP-058: replayed exact $[[186,10,14]]$ promoted class; seven exact classes audited. |
| `results/certificates/exp060_210_18_8_distance.json` | 3,833 B | EXP-060: replayed exact $[[210,18,8]]$ promoted class; thirteen exact classes. |
| `results/certificates/exp064_210_24_4_distance.json` | 3,536 B | EXP-064: replayed exact $[[210,24,4]]$. |
| `results/certificates/exp064_210_14_12_distance.json` | 3,827 B | EXP-064: replayed exact $[[210,14,12]]$. |
| `results/certificates/exp064_210_10_16_distance.json` | 3,966 B | EXP-064: replayed exact $[[210,10,16]]$. |
| `results/processed/exp059_coprime_transport.json` | 1,287 B | Exact all-monomial CRT transport for the three coprime $N=105$ lattices. |
| `results/processed/exp060_n210_ratchet.json` | 54,307 B | Adaptive exact results plus one explicitly witness-only dominated $k=8$ row. |
| `results/certificates/exp067_234_8_18_bundle19_distance.json` | 7,027 B | EXP-067: exact Liang Table III $[[234,8,18]]$ row 2, two rooted cap-16 runs plus replays and physical $X/Z$ witnesses. |
| `results/certificates/exp067_234_8_18_bundle22_distance.json` | 7,024 B | EXP-067: exact Liang Table III $[[234,8,18]]$ row 1, two rooted cap-16 runs plus replays and physical $X/Z$ witnesses. |
| `results/partial_runs/exp067_n234_cluster/` | ~776 KiB | Eight rooted connected-cluster run records, exact MatrixMarket inputs, race-free arm64 solver binary, build manifest, license, and complete dist-m4ri/M4RI source archives. |
| `results/partial_runs/exp068_screen_witnesses/` | ~72 KiB | 57 hash-bound physical fallback witnesses: 29 independent deterministic coset reductions and 28 exact CRT transports, all checked against rebuilt matrices. |
| `results/processed/exp055_odd_lattice_screen.json` | 7,147,344 B | Validator-v11 exact fixed-point screen through $n=234$: 4,862 classes / 150,581 pairs, all 4,658 referenced dominated with record-level proofs, 204 no-reference, zero survivor/undecided. |
| `results/processed/exp066_n234_frontier.json` | 10,398 B | EXP-066 576-map audit refreshed after EXP-067 promotion: 182 hard classes to 30 bundles, all 842 local classes dominated. |
| `results/certificates/exp055_odd_lattice_survivors.json` | 356 B | Empty-survivor certificate hash-bound to the exact $n=234$ screen. |
| `reports/claims_matrix.md` | 35,340 B | Machine-auditable 57-claim provenance ledger: every number → certifying artifact → regenerating command → guarding test (38 VERIFIED / 19 UNGUARDED, status-as-printed in the ledger). |
| `notes/open_status_2026-08-17.md` | 10,729 B | Dated novelty/open-status annex: literature-search audit trail supporting the "to our knowledge new here" list (paper §9). |
| `results/partial_runs/exp039_nogo_module.json` | 66,430 B | EXP-039 assembled result over the 134 parent T-certificates — 61 family-closed parents, 249 rows capped a priori, falsification-gate summary. **NOTE: the task template listed `results/processed/exp039_nogo_module.json`; the artifact actually lives under `results/partial_runs/` (as in paper §8). No `processed/` copy exists; ship the `partial_runs` path.** |
| `results/partial_runs/exp039/` (134 per-parent files) | ~413 KB total | Per-parent SAT/UNSAT certificates with witness vectors, `T_is_exact`, and per-parent wall times (full replay inputs for the gate). |
| `results/processed/exp036_envelope_check.json` | 46,979 B | EXP-036/037: all 7 certified reversals are CSS-dominated at equal n over 162 exact candidates. |
| `results/processed/exp040_saturation_probe.json` | 1,009,661 B | EXP-040 small-lattice probe: 64 parents, 26,898 valid perturbations, 496 strict increases — every one on the theorem law, none above it. |
| `results/processed/exp041_t_crosscheck.json` | 70,288 B | EXP-041: independent T recomputation cross-checking the EXP-039 host values. |

Keep `results/partial_runs/exp039/` in the package only if the total stays under
arXiv's ancillary size budget; the assembled JSON alone already supports the gate.

---

## 5. Reviewer quickstart (regenerate the headline certificates)

Canonical §8 block from the paper (3 commands; **WARNING: its first command is a
10.3-hour full run — see the honest timing ledger**):

```bash
cd math/qec
PYTHONPATH=src .venv/bin/python experiments/exp039_nogo_module.py run --ns 144   # FULL: ~620.7 min sequential
PYTHONPATH=src .venv/bin/python experiments/exp039_nogo_module.py gate
PYTHONPATH=src .venv/bin/python -m pytest tests/test_pbb_nogo.py tests/test_pbb_survival.py -q
```

**Honest timing ledger** (sums of `wall_time_s` over the stored certs
`results/partial_runs/exp039/*.json`, sequential single-thread; a re-run replays the
same SAT searches): n≤36 → 3 certs, <0.01 min; n≤72 → 11 certs, ~0.05 min; n≤108 →
79 certs, ~34.7 min; n≤144 → 117 certs, ~620.7 min; ALL 134 → ~29.0 h. **So paper §8
`run --ns 144` does NOT fit in 30 minutes.** The ≤30-min quickstart therefore substitutes
the smallest honest from-scratch certification pass:

| # | Command (from `math/qec/`) | Expected exit | Runtime | Artifact produced/read |
|---|---|---|---|---|
| 1 | `PYTHONPATH=src OMP_NUM_THREADS=1 .venv/bin/python experiments/exp039_nogo_module.py run --ns 72` | 0 | seconds (11 parents; per-cert wall sum 3 s) | rewrites `results/partial_runs/exp039_nogo_module.json` + per-cert JSONs in `results/partial_runs/exp039/` for the 11 parents n≤72, each with `T_is_exact: true` |
| 2 | `PYTHONPATH=src .venv/bin/python experiments/exp039_nogo_module.py gate` | 0; prints `gate passed: 7/7 certified reversals are permitted`; **1 iff `theorem_refuted: true`** | **1.5 s — measured this session, EXIT=0 measured** | reads/replays stored certs against the 7 certified reversals |
| 3 | `PYTHONPATH=src .venv/bin/python -m pytest tests/test_pbb_theorems.py -q` | 0 (4 passed) | ~1-2 min est. [UNVERIFIED this session — running it needs SAT launches the session CPU policy forbids] | Theorem I legs: Lemmas 2/3 machine-checked on all 202 parents / 368 rows |
| 4 | `PYTHONPATH=src .venv/bin/python -m pytest tests/test_pbb_nogo.py tests/test_pbb_survival.py -q` | 0 (17 + 15 passed) | minutes est. [UNVERIFIED this session — same policy] | Theorem G/H machine checks; replays stored certificates |

Total for steps 1-4: comfortably under 30 min if steps 3-4 run at typical scale;
steps 1-2 alone (fresh SAT-backed recertification of the small parents + falsification
gate) are under 1 minute. Exit-code contract is enforced in the driver:
`exp039_nogo_module.py` `main()` maps `run` to `0` on success, `gate` to `0` on pass /
`1` on any violation, and `assemble` to `1` iff `falsification_gate.theorem_refuted`
(`experiments/exp039_nogo_module.py`, `main()` + `gate_only()`).

Every SAT decision stores its canonical-CNF SHA-256 + encoding version; replay
re-derives verdicts from rebuilt matrices and rejects stamps that fail to hash-bind, so
command 2 can never silently pass a stale artifact (paper §8 mechanism; forgery
regressions in `tests/test_pbb_nogo.py`).

---

## 6. Authors / affiliations / submitter

- **Authors: [TO FILL BY HUMAN]**
- **Affiliations: [TO FILL BY HUMAN]**
- **Contact email: [TO FILL BY HUMAN]**
- **Submitter: [TO FILL BY HUMAN]**

arXiv requires a *named human* arXiv account holder as submitter; an agent/LLM cannot
hold author or submitter status. The compiled `reports/paper_pbb_nogo.tex` currently
ships `\author{}` empty (line 45) — this MUST be filled before any submission. New
submitters may additionally need quant-ph endorsement; budget time for that before the
target date.

---

## 7. License

Choice: **arXiv default license** (non-exclusive license to distribute; grants readers
no redistribution/derivative rights) vs **CC-BY 4.0**.

**Recommendation: CC-BY 4.0.** Rationale: (i) this submission carries substantive
supplements — the claims ledger, JSON certificates, the dated novelty annex — that
readers will want to quote, re-host, and build on; CC-BY grants that reuse
unambiguously with attribution, while the default license leaves every such reuse
legally grey; (ii) the paper's selling point is machine-checkability and open reviewer
access, which a restrictive distribution license weakens; (iii) CC-BY fits downstream
QEC code/data reuse norms and costs a theorem-result paper nothing. If a later journal
venue demands copyright transfer, the arXiv default is the conservative fallback — but
for a no-go result whose value IS the certificate trail, CC-BY maximises adoption.

---

## 8. Pre-submission checklist

Each box names the satisfying test/artifact (paths relative to `math/qec/`).

- [x] **Abstract is byte-verbatim the paper's current abstract, incl. both boxed
  results.** Verbatim-substring check re-run this session: PASS (§9; source:
  `reports/paper_pbb_nogo.md` lines 13-63).
- [x] **Page/table counts in §3 measured from the compiled PDF, not the md.**
  `reports/paper_pbb_nogo.log` → "Output written on paper_pbb_nogo.pdf (18 pages,
  716216 bytes)"; `reports/paper_pbb_nogo.tex` → 7 tables / 0 figures (grep counts).
- [x] **Theorem G machine-checked.** `tests/test_pbb_survival.py` (15 checks) — guard
  for C01-C04 in `reports/claims_matrix.md`.
- [x] **Theorem H + exact T enumeration guarded.** `tests/test_pbb_nogo.py` (17
  checks); artifacts `results/partial_runs/exp039_nogo_module.json` +
  `results/partial_runs/exp039/*.json` (134 certs, all `T_is_exact: true`).
- [x] **Theorem I (forced saturation) legs machine-checked.**
  `tests/test_pbb_theorems.py` (4 checks: Lemma 2 all-202 parents, Lemma 3 all-368
  rows).
- [x] **Falsification gate green on demand.** `exp039_nogo_module.py gate` run this
  session: exit 0, 7/7 certified reversals permitted, slack 0 on all, `theorem_refuted:
  false`.
- [x] **61 family-closed / 249 rows capped / exact at n≤144 claims match the
  artifact.** `results/partial_runs/exp039_nogo_module.json` keys
  `family_closed_parents`, `rows_capped_a_priori`; per-n cert rebucketing (§5 ledger)
  agrees with paper §5 prose (61 closed, 117 certs at n≤144, 134 total).
- [x] **Envelope-domination claim (§6).** `results/processed/exp036_envelope_check.json`
  — all 7 reversals dominated by 162 exact CSS candidates (claims matrix C16/C17).
- [x] **T values independently cross-checked.** `results/processed/exp041_t_crosscheck.json`
  (claims matrix C56/C57).
- [x] **Every quantitative claim lands on a file.** `reports/claims_matrix.md`: 38
  VERIFIED / 19 UNGUARDED of 57, with regenerating commands and guard tests per row.
- [x] **Novelty statement current.** `notes/open_status_2026-08-17.md` dated
  search audit (2026-08-17) supports paper §9 novelty list.
- [x] **Hostile review disposition.** Latest hostile review: NO_CRITICAL_FINDING;
  LOW-1 fixed in the current md. [INFERENCE — disposition conveyed via session
  context; the review report itself is not vendored under `reports/` here.]
- [x] **Sources compile clean.** pdflatex log of 2026-08-17: output written, no
  error lines, `bbl` present (compile artifacts `paper_pbb_nogo.aux/.bbl/.log/.out`).
- [ ] **Author block filled (§6).** Blocks submission: arXiv requires a named human
  submitter.
- [ ] **Human picks license (§7) and submits.** This kit only stages; it posts nothing.

---

## 9. Verification record (abstract verbatim check)

Method (this session): locate the `## Abstract` body in `reports/paper_pbb_nogo.md`,
slice lines 13-63 out byte-exact, embed them above inside a fenced
code block (no quoting/mangling), then verify on the round trip that extracting the
fence contents reproduces the source slice byte-for-byte:

- exact substring of source md: **PASS**
- abstract occurs exactly once in the source md: **PASS**
- embedded copy reconstructs to the source slice byte-for-byte: **PASS**
- boxed-result count in the embedded copy: **2** (`d_Q > d_Z(P) ⇒ k_Q ≤ k_P − T(P)`;
  `d_Q > d_Z(P) ∧ T(P) ≥ k_P/2 ⇒ k_Q = k_P − T(P) exactly`)
