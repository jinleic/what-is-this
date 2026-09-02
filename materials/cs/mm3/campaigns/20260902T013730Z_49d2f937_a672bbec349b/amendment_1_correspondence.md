# Amendment 1 — widen the Route-P product correspondence to the convention space

**Created 2026-09-02, AFTER the first Route-P attempt failed and BEFORE any
amended compute.** Committed path-scoped before the amended run. Campaign:
`20260902T013730Z_49d2f937_a672bbec349b`, gate `mws59-bit-58vs59`,
prereg sha `aede812c…` (unchanged; this amendment only widens one definition
inside Route P and is fully disclosed).

## What happened (preserved verbatim)

The first Route-P run (`verdict_preamendment.json`, `runner_output.json` of
that attempt) reported:

```
P recount of printed Table 2: left=15 right=15 output=29 total=59
P FAILED: no product correspondence between Table 2 and Table 3
G verdict: INCONCLUSIVE-ROUTE-P-FAILED
```

So the **recount succeeded and matched both the paper's claim and the frozen
session-3 machine tally (15/15/29 = 59)**, but the *correspondence* step —
matching Table-2's product numbering to the Table-3 factor blocks — found no
solution inside the space I had pre-registered.

## Why the registered definition was too narrow (defect in my instrument, not a fact about the objects)

pre_statement §3.2 fixed the correspondence space as: a bijection π with
`U2[r] = ε_r·U3[π(r)]`, `V2[r] = ε′_r·V3[π(r)]` and **ε_r·ε′_r = +1**, with the
output map matching under π **without a compensating sign**, and with the A/B/C
coordinate flattenings assumed identical between the two tables. Three of those
constraints are unjustified:

1. **ε·ε′ = −1 is legitimate.** A product may be printed as
   `M_r = (−u)(−v)` or as `(−u)(v)`; in the latter case the product VALUE is
   negated, and the scheme stays valid because the output map's column for that
   product carries the compensating sign. Requiring ε·ε′ = +1 excludes valid
   matchings; the correct requirement is that the output column match with sign
   `ε_r·ε′_r`.
2. **The two tables need not share one 3×3 flattening.** Table 3 is the
   "format readable by the implementation" (the Perminov/SW file format, rows =
   matrix entries, columns = products); Table 2 uses the paper's compact
   `A0..A8`/`B0..B8` indexing. The frozen evidence already documents exactly
   this class of convention gap: gate A had to apply the C-permutation
   `[0,3,6,1,4,7,2,5,8]` to make Perminov's JSON agree with the paper, and
   session 5 recorded a σ² transposition-convention fix. The same involution
   must therefore be admitted here.
3. Coordinate maps may differ **per side** (A, B and the output index are three
   independent flattenings in these file formats).

## The amended correspondence space (finite, tiny, fully verified — no fitting knob)

Admit `σ_A, σ_B, σ_C ∈ {id, T}` with `T = [0,3,6,1,4,7,2,5,8]` (the frozen
transpose involution — no other permutation is admitted), a bijection π on the
23 products, and per-product signs `ε_r, ε′_r ∈ {±1}`. The correspondence is
ACCEPTED only if ALL of the following hold exactly (integer equality, no
tolerance):

- `lrows[r][σ_A(k)] = ε_r · U3[π(r)][k]` for all r, k;
- `rrows[r][σ_B(k)] = ε′_r · V3[π(r)][k]` for all r, k;
- `outs2[k][r] = ε_r·ε′_r · W3[π(r)][σ_C(k)]` for all k, r.

Search order is deterministic (σ_A, σ_B, σ_C in the fixed order
(id,id,id), (id,id,T), …, (T,T,T); products in ascending index; the first
consistent match wins). **This cannot manufacture a false positive:** the
acceptance test is exact equality of every one of the 23·9 + 23·9 + 9·23
entries, and the campaign additionally re-verifies the assembled circuit by
exact expansion and the decomposition by all 729 Brent identities over Z in two
arithmetics. A wrong convention fails on the first row.

## What this amendment does NOT change

The question, the model, the addition-count semantics, the frozen inherited
facts, the budget, the control battery (including the new R3
wrong-correspondence plant, which is *strengthened*: the plant now swaps two
labels inside the widened space and must still be rejected), the verdict rule
and its branch definitions, and the rule that a budget-bound Route-S prefix may
claim only its decided prefix. Route S is unchanged and remains the fallback if
the widened Route P also fails.

## Disclosure

The failed narrow attempt is preserved (`verdict_preamendment.json`,
`runner_output_preamendment.json`, and the log line quoted above). No result
from it is used. The amended Route P runs only after this file is committed.
