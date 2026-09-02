# Clause-by-clause hypothesis check — theorem 4.1 applied to 3-slice
# L-families of T_O; corollary 4.2 applied to (tau control and modules)

Campaign: `2026-09-01T03:21:00Z_routeA_prov`. Written after
theorem_import.md; every check below is executed or reasoned exactly as
labeled. Script: `routeA_clause_check.py` (frozen in this campaign),
output `clause_check.out`.

## Hierarchy of what gets checked

Theorem 4.1 (as printed, with plain-R inequality (4.2) as the
field-independent core) applied to a 3-slice tensor t = (A, B, C):

- H-clause-1 [STRUCTURAL]: tensor has exactly 3 slices with dim U = dim V
  = n, dim W = 3.
- H-clause-2 [ARITHMETIC/STRUCTURAL]: hypothesis "A invertible" — one of
  the three slices must be invertible.
- H-clause-3 [ARITHMETIC]: conclusion consumes the matrix
  B A^{-1} C - C A^{-1} B computed over the base field (here R, exactly);
  the certified bound is then
      rank R(t) ⩾ n + (1/2) rk(BA⁻¹C − CA⁻¹B)
  using proof-step (4.2) (plain rank), with the printed final step
  (border rank via Zariski closure) NOT needed for our real-rank use.

Corollary 4.2 (module form):

- H-clause-4 [STRUCTURAL]: Λ = finite-dimensional associative algebra
  with unity over k, N a unitary left Λ-module, dim N = n; the structural
  tensor of the module structure map is the object.
- H-clause-5 [ARITHMETIC]: rank_N(bc − cb) rank of left-multiplication
  commutator; bound R̲(N) ⩾ n + ½ rank_N(bc − cb).

## APPLICATION TO OUR TENSOR CLASS — the decisive check

Class D2 of the pre-statement: t = (L_u, L_v, L_w), u, v, w linearly
independent in R^8, viewed as an 8 × 8 × 3 tensor (dim U = dim V = 8,
dim W = 3).

Clause-1: t has three 8×8 slices; dim U = dim V = 8 = n; dim W = 3.
SATISFIED (dimension count, exact).

Clause-2: A invertible? For u ≠ 0 in the division algebra O, L_u is
invertible with L_u^{-1} = N(u)^{-1} L_{ū} (owner-verified identities
L_{ū}L_v + L_{v̄}L_u = 2⟨u,v⟩I, L_{ū}L_u = L_u L_{ū} = N(u)I).
SATISFIED exactly (fmpq check on basis triples; structurally for all
nonzero u since N(u) = ||u||² > 0 over R).

Clause-3: the theorem's matrix is M = B A^{-1} C - C A^{-1} B with
A = L_u, B = L_v, C = L_w:
    M = L_v L_u^{-1} L_w - L_w L_u^{-1} L_v.
Route AF's frozen exact identity (commutator_derivation.txt (4)) pins the
RANK of the X,Y-order commutator
    D = [L_{ū}L_v, L_{ū}L_w] = L_{ū}L_vL_{ū}L_w − L_{ū}L_wL_{ū}L_v
via D² = −4 N(u) det Gram(u,v,w) I with det Gram > 0 on independent
triples ⇒ rank D = 8 exactly, universally. Conjugacy pinned exactly:
L_v L_u^{-1} = L_u (L_u^{-1} L_v) L_u^{-1} and likewise for w, so
M = L_u D' L_u^{-1} where D' = [L_u^{-1}L_v, L_u^{-1}L_w]; wait — order:
M = L_v L_u^{-1} L_w − L_w L_u^{-1} L_v = L_u( L_u^{-1}L_v L_u^{-1}L_w −
L_u^{-1}L_w L_u^{-1}L_v )L_u^{-1} = L_u D' L_u^{-1}? Verify:
L_u (L_u^{-1}L_v L_u^{-1}L_w) L_u^{-1} = L_v L_u^{-1} L_w. YES. So
M = L_u D' L_u^{-1}, rank M = rank D'. And D' has the same square
identity: D' = N(u)^{-2} D? Check: L_u^{-1} = N(u)^{-1}L_{ū}, so
L_u^{-1}L_v L_u^{-1}L_w = N(u)^{-2} L_{ū}L_v L_{ū}L_w, hence
D' = N(u)^{-2} D, rank preserved. THEREFORE rank M = 8 exactly on EVERY
independent triple — UNIVERSALLY, by the frozen human proof + exact
coefficient anchors (already adjudicated in Route AF as
MACHINE-VERIFIED identity checks + HUMAN proof deduction).

The theorem then yields, for EVERY linearly independent (u, v, w):

    rank_R(L_u, L_v, L_w) ⩾ 8 + (1/2)(8) = 12.

NOT 13, NOT 16, NOT 14. Exactly 12. This is WEAKER than the already-
certified chain floor 13 (1 + pencil 12, C4-tight). So Route A closes
with the primary text: the best 3-slice commutator bound is 12 < 13 —
the universal commutator-rank-8 theorem (now fully provenance-clean)
feeds a bound that cannot lift S3.

CONTROL check (clause A3, mandatory): tau = (L_1, L_i, L_j) at n = 4:
A = L_1 = I_4 invertible; M = L_i L_j − L_j L_i (with L_1⁻¹ = I);
exact quaternion algebra: ij = k, ji = −k, so L_i L_j − L_j L_i =
L_k − L_{−k} = 2 L_k, rank 4 (L_k invertible). Theorem gives
4 + (1/2)(4) = 6 ⩽ true rank 7. CONSISTENT, non-tight, as Route AF
already recorded. A reading "4 + 4 = 8" would REFUTE the control; no
such reading exists in the primary text (confirming the Route AF
retraction of the naked twin).

Corollary 4.2 applied to (Λ, N) = (H' as algebra, N = H' as module)?
OCTONION CAVEAT (this is where the campaign must be careful): for the
octonion L-family the corollary does NOT apply with Λ = O because O is
NOT associative — the corollary requires "finite dimensional associative
algebra". The way the theorem reaches our non-associative case is
exclusively through the 3-slice TENSOR statement (Theorem 4.1), which
needs no algebra structure on the slices. This is why the S3-family
application must cite Theorem 4.1, not Corollary 4.2. For the tau
CONTROL (Λ = quaternion algebra H', associative, N = H'): corollary
gives R̲(H') ⩾ 4 + ½ rank_{H'}(i·j − j·i) = 4 + ½ rank_H'(2k) = 4 + 2 = 6
for the module N = H' — i.e. R̲_k over algebraically closed k of the
QUATERNION-ALGEBRA module tensor ⩾ 6; the published R_R(tau) = 7 stays
untouched (7 ≥ 6 consistent). [This corollary-4.2 number is recorded for
completeness; NOT load-bearing.]

The specific tau ⊠ s 8×8×3 object (D4, secondary control): slices
blockdiag(τ_p, τ_p); its (A, B, C) = (L_1^(8), L_i^(8), L_j^(8)) with the
8×8 block-diagonal L's — same theorem, M = L_i L_j − L_j L_i =
blockdiag(2L_k, 2L_k), rank 8; bound 8 + 4 = 12 ⩽ 13 ≤ its true rank in
[13,14]. Consistent; nothing new.

## Clause check labels

- Clause-1 (3-slice shape): HUMAN-AUDITED (definitional).
- Clause-2 (A invertible, N(u) > 0): MACHINE-VERIFIED arithmetic on basis
  representatives + owner-verified algebra identities (HUMAN-AUDITED
  deduction for the universal statement).
- Clause-3 (rank M = 8 universal): the deduction is HUMAN-AUDITED
  (commutator_derivation.txt), anchored by MACHINE-VERIFIED exact
  coefficient identity checks (8 adjoint + 64 polar + 56 rank + 28
  singular-plant + quaternion control) — re-run in this campaign
  (clause_check.out) to confirm the frozen artifacts reproduce.
- tau control arithmetic: MACHINE-VERIFIED (exact fmpq, this campaign's
  script).
- Theorem import: CITED-DEPENDENCY (primary read; scans frozen + hashed).
- Consequent: r ⩾ 12 on D2 universally: HUMAN-AUDITED deduction resting
  on CITED-DEPENDENCY theorem + MACHINE-VERIFIED arithmetic of rank M.

## Route A adjudicated outcome

    For all linearly independent (u,v,w) ∈ R^8:
        rank_R(L_u, L_v, L_w) ⩾ 12   (Strassen Thm 4.1/(4.2) with ½)
    — strictly weaker than the standing certified floor 13.
    No form of the Strassen machinery yields 14: the primary text contains
    no 3-slice bound stronger than n + ½ rank(comm) (first-hand negative).

Route A is therefore CLOSED-NEGATIVE for the ≥16 goal at the
provenance layer too: not only does the /2 factor cap the reading at 12,
the naked "n + rank" form does not exist in Strassen 1983. With S3's
well-documented obstruction (13 chain floor +1 gap; Route AF's exact
commutator identity), no live route to 14 from this theorem family. The
16-lower-bound route is dead as a theorem-family route; any revival needs
a DIFFERENT theorem (none found in the primary; Blaeser absent from the
paper; Lickteig not in it).

The published window 18 ≤ R_R(T_O) ≤ 25 is untouched by every number in
this campaign (12 < 13 ≤ the floors feeding 18; no upper-side claim).
Nothing here refutes published work; no escalation needed.
