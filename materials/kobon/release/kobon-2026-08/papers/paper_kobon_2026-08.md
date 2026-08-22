# Capacity, Escape, and a Verified 54: the Kobon triangle problem in late 2026-08

**Working draft, 2026-08-21.** Consolidated main results of the current
computational campaign around the Kobon triangle problem. Companion full
document: `math/kobon/paper_capacity.md`. All empirical statements carry
independent exact verification; all solver-dependant statements are labeled
by their evidence class.

---

## 1. The problem, two conventions

For $n$ distinct affine lines, a **Kobon triangle** is any triangle whose
sides lie on three of the lines (vertices: three pairwise crossing points);
the problem asks for the maximum number of Kobon triangles with pairwise
disjoint open interiors — a **selected family** $S$. Lines not supporting
such a triangle may cross its interior: each such line charges one
crossing incidence to the global count $C=\sum_{T\in S} c_T$, the quantity
that appears on the left of the capacity inequality.

We use the **broad convention** $K_{\rm gen}$: parallels, multipoints
(pencils and coincident crossings), and arbitrary selections of interior
disjoint triangles are all permitted — the convention governed by OEIS
A006066's open rows; our proofs and solvers treat it exactly. The classical
**face convention** $a_3^s$ (Blanc's and BBL's theorem families) counts
triangular faces of simple affine pseudoline arrangements. These agree at
generic-pose extremal examples but differ formally when degeneracy is
exploited; the campaign's lower-bound certificates are stated in $K_{\rm gen}$
and only contextualized against $a_3^s$.

## 2. The capacity theorem (new)

**Theorem 1 ($C$-refined capacity bound).** For any essential arrangement
of $n$ affine lines with $Q$ parallel pairs and multipoints $p$ of
multiplicity $k_p$, and any pairwise-interior-disjoint selection $S$ with
$C$ line-crossing incidences:

$$3|S|+C + 2Q + \sum\nolimits_p k_p(k_p-4) \;\le\; n\,(n-2).$$

- Sharp at $n=9$ (21) and $n=15$ (65); one-off at $n\in\{5,6,7,11,12,13\}$
  of the realized extremal census; correct at every tested instance
  (5.9M machine checks; independent adversarial audit).
- Refines both the unsigned gap bound and the crossing-only forms by exactly
  the two terms necessary and sufficient: pure forms $\sum_p \binom{k_p-2}{2}$
  and $(k_p-2)^2$ are **false** (exhibited exact counter-witnesses); the
  surviving penalty $k_p(k_p-4)$ is exactly tight: $-3$ for a triple point,
  $0$ for $k=4$, $+5$ for $k=5$.
- Proof: per-line injection of side-intervals $\sigma_\ell$ and chords
  $\gamma_\ell$ into $m_\ell+2a_\ell$; the incidence identity (Lemma 7)
  $m_\ell+2a_\ell = n-2-q_\ell-\sum_{p \ni \ell} (k_p-4)$ and the double
  counting $\sum_\ell \sum_{p \ni \ell}(k_p-4)=\sum_p k_p (k_p-4)$ close the
  bound. No degenerate case slips through: the all-parallel class is the
  unique counterexample (Lemma 8).
- Equality regime: at $\Delta=n(n-2)-3T\in\{0,1\}$ the per-line budget is
  essentially saturated, which forces a rigid structural anatomy (every gap
  used; reaching/double-extremeness constraints) — exploited in §4.

## 3. Corollaries: the discriminant view

Rewritten for a fixed target $T$: any $T$-family must satisfy
$C + 2Q + \sum_p k_p(k_p-4) \le \Delta$. At each row of OPEN frontiers:

| $n$ | $T$ | $\Delta$ | Budget shape |
|--:|--:|--:|--|
| 11 | 33 | 0 | faces-only branch; $Q=0$, no concurrency; equality forces exact counts |
| 12 | 39 | 3 | $Q\le1$, or $Q=0$ with $N_3=1$, $C\le3$ |
| 14 | 54 | 6 | $Q\le3$ (no-multipoint branch); $Q\ge4$ needs compensating triple points |
| 18 | 94 | 6 | as $n=14$ |
| 20 | 117 | 9 | $Q\le4$ no-multipoint |
| 14 | **55** | **3** | new broad upper frontier after §5 |

## 4. Exploitability at knapsack-slack zero: Theorem M and its $n=18$ analog

**Theorem M ($n=14$, $Q=3$ branch).** A 14-line arrangement with exactly
three parallel pairs or one parallel triad, no concurrent triples of
lines, cannot carry 54 interior-disjoint triangles.

*Proof architecture.* At $(n,Q,T)=(14,3,54)$ the slack vanishes; the
equality regime forces the *escape/collinear-double-extreme multiplicities*
$m_i \in \{0,1,2\}$ to solve $m_{i-1}+m_i = 2s_i$ on a direction cycle of
class sizes $(s_1,\ldots,s_R)$ with upper and lower bound $0\le m_i \le s_i
s_{i+1}$.

- *Paired layout* [2,2,2]+8 singles: $R=11$ (odd) — the system has a unique
  integer candidate per placement; **all $\binom{11}{3}=165$ placements
  generate a boundary violation** (machine-checked, exact integer solver).
- *Triad layout* [3]+11 singles: $R=12$ (even) — the alternating-sum
  consistency condition equals $\pm4$ at every placement of the 3-class;
  **no rational solution at all** (machine-checked).

*n=18 analog (conditional layers named).* Same mechanism: all
$\binom{15}{3}=455$ paired placements are M1/M2-infeasible (exact
integer solver; gap subtotals 140/90/225, dihedral orbit count 19, Burnside
reflection quotient 231), and all 16 triad placements carry alternating sum
$\pm 4$ (reflection quotient 9). All finite enumerations machine-checked by
`scratch/kobon/n18/appendix_check_n18.py`, solver functions verbatim-copied
from the audited cycle certificate; the geometric lemmas (E0,L1-L4) are
hand-proved in `scratch/kobon/n18/obstruction_n18.md`.

*n=20 analog (conditional).* At $(20,3,117)$ the slack leaves 3 units
($E+C=3$), so the exact equalities no longer follow: a defect-tolerant
extension is the named open lemma (`scratch/kobon/n20/obstruction_n20.md`).

## 5. The headline: $K_{\rm gen}(14)\ge 54$, verified three ways

**Theorem (lower certificate).** $K_{\rm gen}(14)\ge54$.

The witness is Maiorana's exact-rational 14-line arrangement
(`github.com/rufio72/kobon_triangles_k14 @e47c7cfd`, `sol1/`,
SHA256 `83fc26666d…`). This workspace verified:

1. *SHA-pin + declared structure*: file hash matches the repo manifest; 14
   distinct lines; computed concurrency census equals
   `declared_triple_points` $=(0,4,10),(0,11,12)$ (two shared-line triple
   points, $Q=0$).
2. *Two independent exact-rational census implementations* (record
   vertex-separation convention and the stricter cevian-aware open-interior
   test) both count **54**.
3. *The campaign's audited engine* (`math/kobon/engine.py::verify_selection`,
   minimum=54) passes on the exact constructed selection.
4. *All 15 solutions* sol1..sol15 pass the replayable sweep
   (`scratch/kobon/n14/sweep_all_verify.py`, exit 0; persisted inputs +
   manifest in `scratch/kobon/n14/maiorana/`).

Scope discipline: (i) the source README marks the candidates "NOT yet
independently confirmed" — *this workspace is the independent confirmation*;
(ii) OEIS A006066's approved row remains `14 >= 53 54 [Bader]`; the 54 row
is an unapproved OEIS-history proposal — we report local verification and
OEIS state separately; (iii) we do NOT assert the face-convention equality
$K(14)=54$: the bounds in the literature (BBL upper 54) are simple-face
values and the witness is non-simple. The broad upper frontier is now
target $T=55$ ($\Delta=3$, all-degeneracy monolith in flight).

## 6. Machinery (why the certificates are trustworthy)

- **Encoding.** Slope-sorted lines; variables P (cross vs parallel), C
  (multi-point triple), X (line-crossing order), selection S, plus
  triangle-crossing indicators coupled to per-line face bounds; exact totalizers.
- **Verification discipline.** Each UNSAT theorem is only claimed once
  DRAT-verified (`drat-trim -w`); each SAT model becomes a claim only after
  exact-rational realization re-verification (`verify_selection`).
  Sentinel errors found and closed during the campaign: proof-output
  redirection (stdout vs file), per-line identity stray $k_p$ (Lemma 7),
  kernel-state cross-solution contamination in ad-hoc verification (fixed by
  SHA-pinning and single-object fail-closed replay), wrong monolith polarity
  by planting units (regression now asserted).
- **Scale.** $n=10$ prior chain: $K_{\rm gen}(10)=25$ — 11-cube cover, 9
  complete DRAT proofs byte-re-verified, 2 truncated cubes re-proven via 16
  byte-verified subcubes. $n=12,T=39$: 24k-var/2.76M-clause monolith in flight.
  $n=14,T=54$: 1.36M-var/7.84M-clause family plus cube tree. $n=18,T=94$:
  1.6–1.7M-var / 29.8M-clause cubes in flight (archived headers). $n=20,T=117$
  q0 cube: archived header 3,455,871 vars / 59,195,137 clauses.

## 7. Literature pins (verified citations)

- Blanc, J. *The best polynomial bounds for the number of triangles in a
  simple arrangement of $n$ pseudo-lines*. **Geombinatorics 21 (2011),
  no. 1, 5–14.** Preprint arXiv:0801.2845 (2008).
- Bartholdi, N., Blanc, J., Loisel, S. *On simple arrangements… with the
  maximum number of triangles*. Contemp. Math. 453 (2008) 105–116; arXiv:0706.0723.
- Clément, G., Bader, J. *Tighter Upper Bound for the Number of Kobon Triangles*.
  ETH draft, 2007. (Contains the $(n-1)^2/3$ typo; MathWorld fixes it to
  $(n^2-2n-2)/3$.)
- Savchuk, P. *Constructing Optimal Kobon Triangle Arrangements via Table
  Encoding, SAT Solving, and Heuristic Straightening*. arXiv:2507.07951 (2025).
- Bader, J. ETH configuration page; 14-line/53-triangle record.
- OEIS A006066 (current approved rows; table convention `>=` = lower bound).

## 8. Limitations

1. Necessity only: condition (1) bounds and constrains but doesn't close
   exact values by itself; the branch-scoped statements (e.g. $n=11$
   faces-only, $n=14$ $Q\le3$ *in the no-multipoint branch*) are filters, not
   proofs of existence/nonexistence.
2. No mod-6/parity strengthening inside the capacity argument; the $C$
   term captures a different direction.
3. Straight lines only — the machinery is slope/rational intersection, not
   pseudoline stretchability.
4. Deepest open branch: $n=20$ defect-tolerant extension of Theorem M
   ($E+C=3$ slack).
5. Broad upper frontier at $n=14$ is target 55 (open).

*This is a working consolidation; section numbering and exact ledger of
prior $n=9/10$ sealed chains live in `math/kobon/report.md`.*
