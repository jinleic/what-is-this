# Route F status note — written at Route A close (no compute this session)

Route A resolved (closed-negative, see VERDICT.md). Per the assignment,
Route F was next "if Route A resolves or hard-blocks". Route A resolved,
BUT the session compute budget was exhausted on the provenance gate and
the clause-check freeze (budget notice at 200 soft-cap requests).

Honest stance, unchanged owner-verified state:
    rank((L_1, L_i, L_j)) ∈ {13, 14} — OPEN.
    Lower 13: chain 1 + pencil 12 (S3/RouteAF, MACHINE-VERIFIED replays).
    Upper 14: blockwise 7 + 7 tau certificate (MACHINE-VERIFIED modulo
    frozen tau facts).
    Numerical rank-13 candidate searches exist (5 CP + 3 extension seeds,
    5000 evals) with NO witness produced — COMPUTATIONAL-EVIDENCE only.

The named next action is unchanged and stands ready to run in a fresh
campaign (pre-statement first, per discipline):

  af-triple13-existence:
  - Downward (settles 13): exact-rational Krawczyk/interval-Newton on the
    20-variable symmetric-real realification of the rank-13 CP system
    (positivedefiniteness of the 20x20 symmetrized Jacobian is the
    hypothesis to check), seeded from the two retrievable numerical
    candidates (rank13_candidate.npz / rank13_extension_candidate.npz,
    frozen in 2026-08-31T08:02:18Z_routeAF/).
  - Upward (settles 14): exact impossibility argument — no instrument
    identified in any campaign so far; parameterized-symmetry or
    substitution-based approaches not yet exhausted, but no candidate
    theorem found (Route AF routeF_theorem_audit.txt: no source-backed
    additivity/multiplicativity covers tau ⊠ s).

No domain was moved, no search re-run, no result manufactured. OPEN at 13
with no witness produced, exactly as the owner's directive states.
