# PRE-REGISTRATION — Run 1 — gate H-MINLINE — falsification branch

Committed: 2026-09-01 (UTC), by agent `RsPe3dMinline` (lifecycle agent label
`RsPe3dMinline`). This file is written and committed BEFORE any compute of this
run. `campaign.py init --gate H-MINLINE --prereg <this file>` is executed only
after a path-scoped git commit of this file; the manifest will bind
`prereg_sha256` to this commit. Main's process order, adopted as binding:
prereg -> path-scoped commit -> init from the committed file -> copy
byte-identically into the minted run dir -> controls/main compute -> freeze ->
close. No parameter-dependent probe precedes the commit; tool/help inspection
only.

Owner (Main, 2026-09-01, hub + `cs/RESULTS.md`). Role of this run: independent
reproduction of reported H-MINLINE counterexample data and clean falsification
of H-MINLINE as a general claim. A corrected theorem, if proved, is exported to
a SEPARATE run (planned gate `H-MINLINE-DGE3`), never folded into this one.

## 1. H-MINLINE, exactly as previously reported (falsification target)

Source: `cs/RESULTS.md` rs-pe3d section ("New owner finding, sharper than the
theorem's clause (3)"), agent-quoted below [REPRODUCED from the ledger]:

> Exhaustive exact censuses at weight exactly `d` return candidate sets that
> are SET-EQUAL to the minimum lines in argmin-`d` directions, with zero
> non-line dependencies: PN4 `35 = 35` over all `C(140,4) = 15,329,615`
> supports, `(17,(4,4,4))` `48 = 48` (all three directions minimum),
> `(41,(4,4,5))` `40 = 40` (directions 0 and 1 only, direction 2 correctly
> absent at `s_2 = 5 > d`). Status: **MACHINE-VERIFIED** as a finite statement
> about those three instances; **CONJECTURE** ("H-MINLINE") as a general
> claim.

General (informal) form under test, fixed NOW:

(H-MINLINE) For a tensor RS instance $C_i=\Lambda_i\mathrm{RS}(S_i,t_i)$,
$i\in[3]$, $d=\min_i d(C_i)$ with $d(C_i)=s_i-t_i+1$, EVERY
weight-$d$ carrier ($S$ with $\dim(V\cap\mathbb F^S)>0$, $|S|=d$,
$V=\sum_i L_i(C_i)$) is the support of a minimum-weight word of some $C_i$,
$i$-embedded (i.e. an "$i$-minimum line": points agreeing off axis $i$, axis-$i$
coordinates forming a $d$-subset of $S_i$), and conversely every such minimum
line is a carrier.

This run tests H-MINLINE on the reported $t=(1,1,1)$ slice at
$(q,s)=(13,(2,2,4))$, where $d(C_i)=s_i$ and the LINE form is equivalent to
the general form: at $t_i=1$ a minimum-weight word of $C_i$ is the indicator
of a full axis-$i$ line, so every "$i$-minimum line" IS a full line.

Pinned evaluation domains and index convention, fixed now:
- $S_i\le\mathbb F_{13}^\star$ multiplicative subgroups, $|S_i|=s_i$.
  PINNED axes: axis 0 = axis 1 = the order-2 subgroup $\{1,12\}=\{1,-1\}$;
  axis 2 = the order-4 subgroup of squares
  $\{2^{2k}\bmod13:k=0..3\}=\{1,4,3,12\}$, stored SORTED-BY-VALUE as
  $[1,3,4,12]$. Axis lists for INDEXING are exactly these sorted lists
  (axis 0/1: $[1,12]$; axis 2: $[1,3,4,12]$). The point index convention for
  lex position $p=(p_0,p_1,p_2)$ (row-major over axes 0,1,2 with axis 2
  fastest) is $p_i=\mathrm{indexOf}(x_i,\text{axis-}i\text{ list})$ with
  axis lists sorted ascending. Both this convention and any alternative
  consistent convention give the same partition into fibers; only the NAMES
  of index tuples differ. Reported index pair `(0,12)` refers to the PRIOR
  census's lexicographic convention (points at flat positions 0 and 12 of
  the 16-point grid); this run reproduces the pair by convention-mapping,
  not by re-adopting the old list order.
- $t_i=1$: $\mathrm{RS}(S_i,1)$ = constants ($\deg f<1$), so $L_i(C_i)$ = line
  space of axis-$i$ constants.

## 2. Instance, recount targets, methods, controls (all fixed now)

Instance E0: $q=13$, $s=(2,2,4)$, $t=(1,1,1)$, $\Lambda=\mathrm{Id}$,
$d=\min_i s_i=2$, $N=16$.

RECOUNT TARGETS (fixed now, from the reported census; exact integers, zero
tolerance — all three must match or the run halts for a bug, not a
refutation):
1. weight-$d=2$ carriers: **24**;
2. among them, axis-line supports: **16** — pinned reconciliation: at
   $t_i=1$ minimum words of $C_i$ have support a FULL axis-$i$ line
   ($i\in\{0,1\}$ at this instance: $s_1s_2\binom{s_0}{2}=8\cdot1=8$ for
   axis 0, $8\cdot1=8$ for axis 1); axis 2 contributes none since a
   4-point axis-2 line has weight $4>d=2$. The prior notation
   "minimum lines in argmin-$d$ directions" therefore coincides with
   {full axis-0 and axis-1 lines} here, and 16 = 8+8 is the line count.
3. off-line diagonal pairs: **8**, including — under the prior census's
   lexicographic flat-position convention — the pair of points with flat
   indices **(0,12)**; this run independently exhibits that exact pair with
   an explicit codeword.

### 2.1 Methods

All arithmetic exact $\mathbb F_{13}$ (Python ints mod 13; numpy int64 for
index arithmetic only; ZERO floating point anywhere).

(a) KERNEL ROUTE. $V=\ker(H_0\otimes H_1\otimes H_2)$ where $H_i$ is a GRS
parity-check for $\mathrm{RS}(S_i,1)$: a single constancy row
$(1,\dots,1)$ on each axis-1D fiber. Membership of a support $S$:
$\dim(V\cap\mathbb F^S)=|S|-\mathrm{rank}(H[:,S])$, $H$ the stacked
Kronecker-check rows restricted to $S$-coordinates (exact mod-13 Gaussian
elimination). Controls compare the kernel route against an INDEPENDENT
REPRESENTATION of $V$ (basis route below).

The instrument re-derives from data: the 8-diagonal classification, and
Main's decomposition $24 = 4\cdot\binom{4}{2}$ — pinned reading: each of the
24 carrier pairs lies entirely inside one axis-2 fiber plate $x_2=c$
($c\in S_2$, 4 plates, each plate a $2\times2$ grid with
$\binom{4}{2}=6$ pairs; the plate partition of the 24 is recorded, with 16
of the 24 full axis-0/axis-1 lines and the other 8 diagonals).


(b) BASIS ROUTE (independent): $V=\sum_i L_i$ computed as a column basis by
exact row reduction of the stack of the three lifted line-class bases, and
$V\cap\mathbb F^S$ tested by solving the coefficient system for the
17 vectors $\{\mathbb 1_S\}\cup\mathrm{basis}(V)$. A support is a RECORDED
carrier only if BOTH routes agree its intersection is nonzero.


(c) CLASSIFICATION of every carrier: LINE (the two points agree in one full
axis: at this instance exactly the full axis-0 or axis-1 lines) or OFF-LINE
diagonal (no axis fully agreed). Supports enumerated exactly: all
$\binom{16}{2}=120$ pairs (weight $d=2$ exactly), full classification plus
counts.

(d) WITNESS DELIVERABLE (acceptance): emit at least one explicit off-line
carrier — a pair agreeing in NO axis — with the exact defining codeword
given as (i) the 16 field values (both pair entries nonzero, all others
zero), (ii) an exact decomposition $v=x_0\otimes a_1\otimes a_2+
x_0'\otimes a_1'\otimes a_2'$-style axis-constant tensor sum AND as a
single tangent vector $g\otimes\mathbb 1\otimes\delta+\mathbb 1\otimes
h\otimes\delta$, and (iii) a kernel-consistency certificate (exact
$H v^{\!\top}=0$). Pre-registered expectation from Main's derivation, TO BE
VERIFIED NOT ASSUMED: with $g=(1,0)$, $h=(0,-1)$ on the plane $x_2=1$ the
tangent vector $v=s_{e_0}\otimes\mathbb 1\otimes\delta_1+\mathbb
1\otimes s_{e_1}\otimes\delta_1$ has $2\times2$ plane matrix
$[[1,0],[0,-1]]$, support the ANTI-diagonal pair — the run computes and
records the actual supported pair under the pinned convention and maps it
to the prior census's flat-(0,12) form.

(e) CONTROLS (all pre-registered, exact pass/fail):
- C1 unique-min: $(13,(2,3,4),(1,1,1))$, $d=2$: every weight-2 carrier is a
  full line; fixed targets **12 line carriers / 0 off-line**
  ($12=s_1s_2\binom{s_0}{2}=4\cdot3\cdot1$, direction 0 argmin-unique).
- C2 tied d=2: $(13,(3,3,4),(1,1,1))$: line carriers in argmin directions =
  $24$, and total weight-2 carriers $=24$ per the report — targets **24/24**;
  any off-line carrier here is a fresh falsification datum (none reported).
- C3 planted noncarrier REJECT: the support test MUST return dim 0 on a
  planted never-carrier — pinned plant at E0: the pair $\{(0,0,0),(1,1,1)\}$
  (differing in all three axes) AND the pair $\{(0,0,0),(0,0,1)\}$
  (half of an axis-2 line, weight below the line's minimum) — both asserted
  dim 0 by the kernel route; the instrument must ALSO accept a known-true
  carrier (any full axis-0 line at E0, dim $\ge1$), and the size-2 extractor
  fed a synthetic vector carried by a known-dependent pair MUST return
  dim $>0$. Known-true accept + nearby-plant reject must both pass before
  E0 data are trusted.

- C4 cross-check: every recorded carrier count must agree between kernel and
  basis routes; any disagreement halts the run (bug, not data).
- C5 fail-loud: any count mismatch vs targets 24/16/8 (E0), 12 (C1), 24 (C2)
  halts falsification output and converts the row to an instrument-defect
  record.

## 4. Adjudication, fixed now

- VERDICT ARM: E0 recount yields 24 carriers of which 16 are lines and 8 are
  off-line diagonals, the explicit witness pair contains flat indices (0,12)
  under the mapped convention, and controls C1-C5 all pass
  => **H-MINLINE FALSIFIED on this instance** (weight-2 carriers exist that
  are not minimum-line supports); terminal verdict for Run 1:
  FROZEN-NEGATIVE.
- Any count mismatch after defect resolution: instrument-defect record; no
  falsification claim.
- If, unexpectedly, E0 returns exactly 16 lines and 0 diagonals: the run
  records the anomaly, halts, and re-examines which slice the report meant
  BEFORE any verdict. No verdict without reproduction of the reported
  counts.
- Labels: every emitted carrier is MACHINE-VERIFIED (exact F_13, both routes);
  the falsification claim is a finite statement about E0/C1/C2 only.

## 5. Proof-branch export (Run 2, separate gate — planned, not decided here)

If a corrected general theorem is proved during this run's analysis window
(circuit characterization; parallel-class hypothesis from Main's steering:
sparks $d_i=s_i-t_i+1$, only factors with at most one spark-2 factor admit the
fiber/circuit characterization; $d_i\ge3$ is the clean sufficient regime),
it is exported to a SEPARATE future campaign (gate `H-MINLINE-DGE3`) with its
own prereg, commit, init, exhaustive controls, freeze, and its own verdict —
per Main's verdict-semantics instruction. This file's verdict applies to
H-MINLINE's line characterization only. Run 2 MUST NOT be opened until Run 1
is frozen and closed.

## 6. Rule-7 scope sentence (fixed now)

Swept: exact $\mathbb F_{13}$ weight-2 census over all $\binom{16}{2}=120$
pairs at E0 $(13,(2,2,4),t=(1,1,1),\Lambda=\mathrm{Id})$, plus C1 $(13,(2,3,4))$
$\binom{24}{2}=276$ pairs, C2 $(13,(3,3,4))$ $\binom{36}{2}=630$ pairs — all
three censuses complete (no sampling anywhere), dual-route (kernel and basis)
agreement mandatory per carrier class.
Not swept, plainly: supports of weight $\ge3$ at any instance; any other $q$,
$s$, $t$, $\eta$; non-identity $\Lambda$; asymptotic claims; Conjecture 4.2
itself; the general circuit theorem (exported to Run 2's domain, not
duplicated here); global $\rho_{\rm inst}$ at any instance.
