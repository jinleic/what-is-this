# Campaign manifest — 2026-08-30T174243Z_edbdd408-840f-42b4-9be0-192a95783ab9

Agent `Mm3RecordAttack`, 2026-08-30. Gate C continuation: the record attack over the
ternary-unimodular (non-monomial) orientation space, following the frozen gate-C sweep
`2026-08-30T113838Z_*`.

## Outcome first

**NO 54-addition scheme exists over the pre-registered named set S; the certified
minimum total addition count over S is exactly 55.** The published record 55
(arXiv:2607.28676) is not beaten by any orientation in S, and every orientation in S
now carries a checker-verified certificate. Nothing was escalated: no record-breaking
scheme was found, and no published claim is contradicted.

## The named set S (pre-registered, exhausted exactly)

- Decompositions: `paper55` (arXiv:2607.28676 §5.2 blocks) and `sun56`
  (arXiv:2604.27645 shipped verify.py) — the two families whose certified total-lb was 55.
- G over all **6960** ternary unimodular 3×3 matrices (|det| = ±1; counted by exact
  brute force over all 3^9 = 19,683 ternary matrices: 11,808 invertible over Z,
  of which 6960 have |det| = 1; 48 of those are the monomial signed permutations).
- Action: the frozen verified isotropy action with **(X, Y, Z) = (G, G, G)**
  (the diagonal slice), composed with the σ-orbit {Id, σ, σ²}.
- Total triples in S: 2 × 6960 × 3 = **41,760**; every one was decided (no sampling,
  no prefix): S is **exhausted exactly** (two independent paths, below).

## Census result (MACHINE-VERIFIED, exact integer + fmpz arithmetic)

Per decomposition, over all 6960 G:

| G class | count | verdict |
|---|---|---|
| non-monomial, leaves ternary alphabet under (G,G,G) | 6912 | factors of the sandwiched triple exit {−1,0,1}^9 — NOT in the ternary-alphabet model; excluded with reason, not silently dropped |
| ternary + Brent-valid | 48 | exactly the monomial signed permutations |

- Zero G kept ternarity but failed Brent (0 in the `ternary_brent_fail` cell of
  both paths' verdicts) — the tensor automorphism never broke validity; the wall
  is purely the alphabet boundary.
- Ternarity is σ-invariant (σ permutes/transposes ternary blocks), so the census is
  σ-complete by structure; the 48 survivors × 3 σ were additionally Brent-checked
  per class (288 checks, 0 failures).
- For the 288 valid orientations: d-counters + complete floor-DFS + transposition
  ⇒ certified total LB = **55 for every class**, reproducing the frozen gate-C
  landscape (min 55). No orientation of S reaches 54.

## Independent verification paths (rule 14 discipline)

1. **Kernel path** (`census_kernel_path.json`): pure-Python integer arithmetic,
   my own re-implementation of sandwich + 729-cell Brent + bound layer.
2. **fmpz path** (`census_fmpz_path.py` → `census_fmpz_path.json`, run under
   `cs/.venv`, python-flint 0.9.0): the **frozen** gate-C sandwich implementation
   (`gatec_sweep.sandwich`) + Brent in fmpz + the frozen bound layer
   (`gate_b_floor.prep`/`subset_dfs`). Both paths agree on
   6912 / 0 / 48 verdicts and on the 55 totals.
3. **Anchors** (frozen harness `verify_anchors.py`): 55/58/56/59/60 all reproduced
   before any new number was trusted.
4. **Frozen-code cross-check**: frozen `d_count` matches my counter on the base
   decompositions (d = 12/13/13 for paper55).

## Checker-verified impossibility certificates (the ≥ 55 upgrade)

All 18 side-instances of the valid orientations ({paper55, sun56} × σ ∈ {0,1,2} ×
{L, R, Ofac}) were decided at T = d(F):

- **15 UNSAT** floor instances → kissat 4.0.4 binary DRAT proofs, converted to
  LRAT by `drat-trim -L`, then **verified by BOTH pinned checkers**
  (`drat-trim` and `lrat-check`, both built from the frozen source snapshot
  `2026-08-30T031544Z_*/tools_snapshot/`, git 2e3b2dc): `s VERIFIED` /
  `c VERIFIED`, exit 0, on every instance. Verbatim output in
  `checker_output/verbatim.txt`. Cross-solver: CaDiCaL 3.0.1 also UNSAT with an
  independently checker-accepted DRAT proof on the audited instance
  (`paper55_s0_L_d12*`).
- **3 floor-achievable** instances (sun56 W-factor d=16 at σ⁰; sun56 R d=16 at σ¹;
  sun56 L d=16 at σ²) — these are theknown achievable Sun floors, matching the
  frozen campaign; they are witnesses, not UNSAT instances, so no proof log exists
  for them; their exactness is witnessed by the frozen 16-gate schedule.
- Certificates re-derived totals: 55 for all six (D, σ) classes — in exact agreement
  with the d-counter/DFS layer and with the frozen certified landscape.

## The two reductions, stated exactly (assignment point 3)

1. **Inside the ternary alphabet**: the ternary-unimodular set contributes only
   **48/6960** usable G per decomposition — reduction factor **6960/48 = 145**
   exactly. This is a census result over the DIAGONAL slice, machine-verified;
   it is NOT the monomial-transfer theorem.
2. **The monomial-transfer boundary**: the frozen theorem covers exactly the 48
   signed-permutation matrices (cost-invariant per orientation). The 6912
   non-monomial matrices are **excluded by the alphabet, not by the theorem** —
   the honest split is: theorem covers 48; census excludes 6912 empirically
   (alphabet exit, instance-exact); the off-diagonal 6960³ × 3 = 1,011,460,608,000
   product remains UNSWEPT. The gate-C agent's framing is inherited verbatim: the
   non-monomial space is a wall, not a sweep.

## Structural finding recorded (rule 5/17 honesty)

The README phrase "the 3×3 ternary unimodular group" is a **misnomer**: the 6960
ternary |detEndpoints|=1 matrices do NOT close under multiplication (counterexamples
in `census_summary.json`: products of unimodular ternary pairs have entries ±2; the
48 monomials do close — hyperoctahedral B₃). The per-element orientation action is
still well-defined (each G is a tensor automorphism with integral inverse), so the
census is sound; only the "group"/orbit language needed correcting. Correction
recorded here rather than silently patched (rule 5).

## Encoding soundness (assignment point 4)

- **No new SAT/ILP encoding was introduced** by this campaign: every decision came
  from (a) the exact d-counter + complete subset-DFS (no solver, no float) and
  (b) the frozen `build_floor_cnf` encoding already audited and checker-certified
  in campaigns `010534Z`/`031544Z` (FORBID-only clauses beyond the gate-B template;
  both-direction model ⇔ circuit proof; no positive-unit under-constraint of the
  arXiv:2607.29291 shape). No new UNSAT rests on a new encoding.
- **HiGHS was NOT used** — no ILP ran, so the 1e-7 tolerance trap is not exercised
  in this campaign. All arithmetic is exact (Python int / fmpz).
- **Ring / F_2 caveat** (assignment point 4): confirmed — an F₂-valid scheme does
  NOT need to satisfy 729/729 over Z, so F₂ schemes are outside this model; S
  contains none (all factors are ternary over Z and every claimed orientation was
  Brent-verified 729/729 over Z in fmpz). arXiv:2607.29291's under-constraining
  defect shape is structurally inexpressible in the floor CNF (no positive unit
  clauses exist beyond the cover clauses that are the theorem's consequences; audit
  inherited from `gate_b_encoding_audit.md`).

## What was NOT searched (mandatory rule-7 scope sentence)

This campaign swept exactly: the two decompositions `paper55` and `sun56`, each
under the full 6960-element ternary-unimodular diagonal slice (X, Y, Z) = (G, G, G)
composed with the 3-element σ-orbit — 41,760 (D, G, σ) triples, each decided
exactly — and it produced checker-verified floor-impossibility certificates for all
15 non-achievable side instances; it did NOT sweep the off-diagonal product
6960³ × 3 = 1,011,460,608,000 (X, Y, Z) triples (pre-registered as budget-out in
`pre_statement.md` §1 and §5), did NOT sweep any decomposition outside the five
named public ones (in particular Laderman, Smirnov's other schemes,
Schwartz–Vaknin, or any unpublished or F₂-only scheme, whose Brent-failure over Z
makes them ineligible for this model), did NOT sweep GL(3,ℚ)/GL(3,ℤ) sandwiches
with non-ternary entries, did NOT run any upper-bound synthesis beyond the known
witnesses, and therefore establishes NO universal no-54 claim: a 54-addition
scheme could still exist outside S (off-diagonal non-monomial sandwiches, other
decompositions, or beyond the ternary alphabet).

## Files

- `pre_statement.md` — committed (git `fb731d4`) before any compute.
- `census_kernel_path.json` / `census_fmpz_path.py` / `census_fmpz_path.json` —
  the two independent census paths.
- `census_summary.json` — verdicts, reduction factors, closure counterexamples.
- `cert_gen.py` — certificate driver (frozen as run, after one repaired syntax
  error and one repaired loop deletion, both pre-execution, recorded here).
- `certificates_result.json` — 18 instances with checker verdicts + SHA-256.
- `cnf/` (15 DIMACS), `proofs/` (15 kissat DRAT + 15 LRAT + CaDiCaL cross-proof),
  `checker_output/verbatim.txt` — solver + both checkers, verbatim, exit codes.
- `checksums.txt`, `versions.txt`.

## Tool versions

- python 3.14.3 (`cs/.venv`), python-flint 0.9.0, sympy 1.14.0.
- kissat 4.0.4, CaDiCaL 3.0.1 (PATH), drat-trim + lrat-check built from
  `tools_snapshot/drat-trim.c` ("Last edit: April 21, 2024", git 2e3b2dc):
  binary SHA-256 drat-trim `111b0405566d55629f5d391b80b3220cc7a3ebc3cd406a35895ff65cc9acd5e4`,
  lrat-check `b4bdebfcc40da664be2fd416451b43f9e764e5e92383ce2b4e9d1f05148e9b07`.
- Machine: macOS 26.5.2, Apple M3 Ultra; `nice -n 10`, single-threaded.
