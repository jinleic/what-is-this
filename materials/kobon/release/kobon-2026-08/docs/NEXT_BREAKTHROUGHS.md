# Next breakthroughs: ranked program (2026-08-21)

Framing correction that drives this list: the Kobon problem's *broad*
convention (Clément–Bader: any triangle bounded by three lines, pairwise
non-overlapping interiors, crossings and degeneracies allowed) is **not** the
face convention in which the sharp asymptotic bounds of BBL/Blanc are proved.
Every published *upper* bound must therefore be classified before use. Our
capacity theorem gives per-line accounting valid in the broad convention;
whether the Clément–Bader mod-6 bound is itself broad is being pinned
verbatim from their draft (agent `NoveltyCheck`) — it decides whether the
open windows below are 1 unit wide or ~13 units wide.

OEIS A006066 was updated **2026-08-21 19:04** to `a(14) = 54 [Maiorana]`,
i.e. 54 is now recorded as *exact*. We independently verified the
construction (three implementations, 15/15 witnesses). The upper half of
that claim is exactly what our in-flight $T=55$ run decides.

---
## 2026-08-24 endpoint-front update

The \(n=12,T=39\) front now has two stronger sound discovery lanes:
`n12-gap-endpoint-mi-k4-t39` (endpoint closure, two required
multipoint-incident faces, direct \(K(4)=2\), and 11-/10-line heredity) and
`n12-gap-endpoint-sub9-t39` (also every checked \(K(9)=21\) cut). They are
375,037 / 1,005,366 and 842,317 / 1,917,926 variables/clauses respectively.
Both are live without proof logging.

The next action is conditional and exact:

1. SAT: decode, straighten, and verify over \(\mathbb Q\); a verified
   39-face arrangement is the breakthrough.
2. UNSAT: rerun the byte-pinned formula with proof logging and independently
   check the proof before claiming \(K(12)=38\).
3. No verdict: cube on the forced multipoint-incident face pair and use the
   endpoint-closure guards as the split surface; do not add the negative
   exact-face or Grassmann--Plücker experiments to the production lane.

---


## Tier 0 — decisive experiments already running

1. **$n=14$, $T=55$ all-degeneracy monolith** (`n14-t55-monolith`).
   *UNSAT* ⇒ with the verified 54 witness this is the first rigorous
   broad-convention proof of the brand-new OEIS entry $a(14)=54$ (promote to
   theorem with DRAT). *SAT* ⇒ **new world record 55** and a refutation of
   the freshly approved entry. Highest value per CPU-hour in the program.
   Follow-up: `n14-t56-monolith` ($\Delta=0$ rung) is queued behind it.
2. **$n=9$, $T=22$ DRAT verification** (`drat-n9-t22`, running) ⇒ certified
   $K_{\rm gen}(9)=21$. Companion: $K_{\rm gen}(8)=15$ is **already
   DRAT-certified** (`s VERIFIED`, 2026-08-21) — the first machine-checkable
   certificate for that OEIS row.
3. **$n=13$, $T=48$ ladder rung** (`ladder-n13-t48`, running). Independent
   falsification test of the Clément–Bader draft bound at $n\equiv1\pmod 6$:
   capacity says a 48-family needs $N_3\ge1$; C–B (if broad) says it cannot
   exist. A SAT model would refute an unpublished-but-cited bound.

## Tier 1 — cheap certified closures (hours–days each)

4. Ladder rungs for every "already exact" small row, in the broad
   convention with all degeneracies free: $n=12/T=39$ and $n=11/T=33$
   (both in flight), then $n=16/T=73$, $n=17/T=86$, $n=19/T=108$,
   $n=15/T=66$, $n=21/T=134$. Each rung is a self-contained
   "is the published value actually exact when crossings are allowed?"
   experiment. Small $n$ are minutes-to-hours (the $n=8$ rung took 75 s and
   its proof verified in 63 s).
5. Promote each discovery UNSAT to a theorem via kissat+`drat-trim`; the
   $n=8$ pipeline (50 MB proof, 63 s verification) is the template and the
   costs scale gently for $n\le13$.

## Tier 2 — the real mathematical prize

6. **Kill the triple-point credit.** Our per-line inequality allows
   $2a_\ell$ tokens, which is what lets $\sum_p k_p(k_p-4)$ go negative and
   opens the ceiling above Tamura. Conjecture: the tokens are never all
   realisable, so the true inequality has penalty $\ge0$ for every
   multiplicity, giving $K_{\rm gen}(n)\le n(n-2)/3$ *with degeneracies*.
   Evidence: in the whole sharpness census, equality occurs **only** at
   $N_3=0$ (attained at $n=3,5,9,15,17$), and every known witness's exact
   optimum equals its published count. Attack: characterise when a triangle
   can corner at a multipoint from both sides simultaneously; the
   double-extreme analysis of Theorem M is the right local tool.
7. **Defect-tolerant obstruction lemma.** `proofs/defect_lemma.md` proves the
   $m$-vector system alone cannot close $n=20$, $T=117$ (explicit feasible
   defect families: paired $b=2$ with runs $(6,0,8)$; triad $b=4$). The
   missing ingredient is an endpoint/location constraint on where a defect
   may sit. Success would close $n=18$ and $n=20$ branches wholesale.
8. **Broad mod-6 strengthening.** Derive the $-1$ correction for
   $n\equiv0,2\pmod 6$ inside our per-line framework (rather than importing
   it), which would make the campaign self-contained and settle the scope
   question permanently.
9. **Lower-order asymptotics.** The slack table
   $\Delta(n)=n(n-2)-3K(n)$ for $n\le19$ reads
   $0,2,0,3,2,3,0,5,3,6,2,6,0,0,2$. Is $\Delta(n)=\Theta(n)$? A construction
   theory for $\Delta$ would answer "how far below Tamura the truth lies",
   the central quantitative question of the problem.

## Tier 3 — constructions (records)

10. **Degeneracy-boosted records at $n=18,20$.** Maiorana's step from 53 to
    54 used two shared-line triple points. The same pattern is untested at
    $n=18$ (published 93; broad window open above it) and $n=20$ (116). Plant
    $N_3=2$ shared-line triples and run target 94/95 and 117/118.
11. **Exact-optimum-guided search.** `verification/max_packing.py` computes
    an arrangement's exact broad optimum in seconds ($n=14$: 6 s; $n=20$:
    30 s). Use it as the objective of a search over exact-rational
    arrangements (neighbourhoods of the 15 verified 54-witnesses, of Bader's
    18/20-line records, and of planted-degeneracy families). This makes
    heuristic search rigorous: every candidate is scored by a certified
    optimum, never by a float count.
12. **Reduction-lemma harvest.** At Tamura-tight $n$ ($n=15,17,19,21,23,27$
    with perfect arrangements known) a triple-point-free maximum family must
    consist of faces, so the published face optima transfer to the broad
    convention; the only gap is arrangements *with* triple points. Combining
    the reduction lemma with item 6 would close all of those rows at once.

## Method notes worth institutionalising

- Every claim carries an evidence class: PROVED / DRAT-CERTIFIED /
  MACHINE-CHECKED-ENUMERATION / DISCOVERY-ONLY / PENDING.
- The exact-rational verifier is the gate: no SAT model is a construction
  until `engine.verify_selection` accepts it, and no UNSAT is a theorem until
  `drat-trim` prints `s VERIFIED`.
- Prefer instances whose *degeneracy pattern is planted*: they are orders of
  magnitude easier than free-degeneracy monoliths and they cover the same
  space when the planting is a proven case split.
