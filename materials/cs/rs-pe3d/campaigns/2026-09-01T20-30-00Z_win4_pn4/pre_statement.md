# Pre-statement — campaign `2026-09-01T20-30-00Z_win4_pn4` (Campaign W, support
# window at PN4)

Committed as the FIRST file of the campaign, before the runner source, the
anchors, the kill-restart demonstration, or any census compute.
Agent: RsPe3dMW. UTC at writing: 2026-09-01T19:58Z.

## 1. Domain (fixed before any compute; Main-ruled option (a), 2026-09-01)

- **Instance:** the frozen P2 row **PN4 = (q=421, s=(4,5,7), t=(1,1,1),
  Λ=(1,1,1) per axis)** — a pairwise-coprime triple (4·5·7 = 140 = N; q−1 =
  420). This is exactly the instance class the P2 battery marked
  UNTESTABLE-in-window (its weight-≤3 window is empty — now a theorem: M's
  Theorem T, first nonempty weight d = min_i s_i = 4), and the instance the
  Assignment named ("the PN4 N=140 class"). Main's dispatch said
  "min(s)≤3 instance" — an owner phrasing error, corrected by Main's
  follow-up ruling: option (a) = PN4 ITSELF at N=140, w=4 slice. The
  alternative option (b) (a hypothetical min(s)≤3 sibling at N=140) is
  NOT run here; it is recorded below as a named follow-up candidate per
  Main's instruction.
- **Window:** point supports of weight EXACTLY 4 (the first nonempty weight:
  |S| = 4, i.e. all C(140,4) = 15,329,615 point 4-subsets). NOT w ≤ 4:
  weights ≤ 3 are known-empty by Theorem T (M campaign, 16/16 controls,
  C-PN4 at this very instance: 457,450/457,450 empty) and add zero
  information — they are excluded to keep the census exact-scope. The
  support weight ≥ 5 window stays OPEN (not run here).
- **Everything else frozen identical to the P2 battery convention:**
  V = span of the three direction line classes (t = (1,1,1): all-ones";
  candidates = normalized reduced basis of each nonzero V ∩ F^S;
  per-witness delta computed EXACTLY by the exhaustive ascending-cost
  line-tuple sweep (all cheaper line-tuple triples proved absent, ≥1 present
  at the claimed cost); ρ^window = min over supports S of wt(S)/delta(S)
  with wt(S) = 4. The census counts every support with V ∩ F^S ≠ {0} and
  records EVERY nonzero-support dimension along with the basis vectors for
  the subsequent per-support delta sweep-adjudication.
- **Adjudication object = rho^window** (best ratio in the w = 4 window),
  exactly as labeled in P2/gate-B — explicitly NOT ρ_inst (the global value
  stays OPEN; heavier supports could in principle beat it; the certified
  interval this campaign will produce is
  ρ_inst ∈ [1/(3N), ρ^window] = [1/420, ρ^window] when ρ^window lands).
- **Predictions (pre-registered):** P1. ρ^window ≤ 1 (a line word gives
  wt 4 = delta ≥ 4 trivially... to be proven, not assumed: the line-word
  delta at PN4 is exactly 4 ONLY IF no cheaper line-tuple closes that word —
  E2' at M verified only its membership in V). P2. Every support S with
  V ∩ F^S ≠ {0} either contains a full axis-0 line (4 points) OR is exactly
  the support of a non-line word (to be tested by the census, not predicted
  strongly). P3 (the P-relevant boundary): under pattern P ("{2,3} ⊆ s ⟺
  ρ^window = 3/5" at the ≤3 window), the PN4 row had NO ≤3-window ratio;
  at the w=4 window P makes NO committed prediction — rho^window_w4 = 3/5
  would be REMARKABLE (P's scope has been ≤3-window-only so far) and is
  flagged as the single most valuable possible outcome; rho^window_w4 = 1
  (all line words) would be the expected floor per the ≤3-window history
  at non-{2,3} rows (PN1-PN9 landed 1 or 1/2, never sub-1/2).
  NO other numeric prediction is committed — this campaign MEASURES and
  reports, with misses prominent if any expectation above fails.

## 2. Instrument: vectorized census kernel + validation ladder

The census kernel is vectorized (numpy): one-time V-basis setup (rref over
all 421 lines' classes: 70+56+60 = 186?? — the actual count is N/s_i lines
per direction: 421/4=105.25 WRONG — recompute: N=140; lines per direction =
N/s_i = 35/28/20, total 83 line classes before rref-dedup; V-basis = 68
rows after rrep). Per support S (encoded as a 140-bit mask), V ∩ F^S has dim
> 0 iff the left-kernel of B_off (V-basis' off-S columns) is nonzero; the
kernel is computed by exact modular-rref of the (N−|S|)=136-column submatrix
transpose. The module uses python-int residues mod q throughout (no float
anywhere); numpy int64 carries residues (q = 421 « 2^31, products ≤ 421²
« 2^63 — no overflow even before mod).

**Validation ladder (MUST pass, in order, before the long run):**
1. **Anchor (rule 17b, startup):** the kernel reproduces the 13/13 frozen
   gate-B tallies (B01..B13 nonzero-support counts, from
   `campaigns/2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB/payload.json`
   `candidate_count_total` where the row's window is the ≤3 census...) —
   precisely: for each of the 13 rows the kernel's w ≤ 3 nonzero-support
   count EQUALS the row's frozen nonzero-support census figure. ALSO the
   two P2 vectorized anchors: PN6 → 3486, PP2 → 917 (same numbers M
   independently re-landed). A mismatch on ANY row halts the launch.
2. **Cross-check vs the unvectorized path at small N:** the full
   `_profile_census` (frozen gate-B path, exact `DeltaEngine rref` per
   support, no numpy) counts nonempty supports at (13,(2,2,4)) w ≤ 3 and
   (5,(2,2,2)) w ≤ 3 (?–use the unit instance) must AGREE with the
   vectorized kernel exactly.
3. **The MN supporting-check at the target instance:** the kernel's w ≤ 3
   per-support dims at PN4 are all zero AND its w = 4 count on a sampled
   subset (every 1,000,000th support in lex order, plus all supports
   sharing ≥3 points with a fixed direction-0 line: 136 + 3·C(136,?) —
   concretely: all C(137,1)+... NO — keep it simple and exact: the w=4
   supports CONTAINING the fixed direction-0 line {L0∪(1 extra point)}
   must be nonempty: 136 of them, one per extra point, line-word +
   extra-cancel... careful: adding a 5th point to a full support of a
   4-point line can force the word to vanish... the DIRECT check: every
   w=4 support that IS a direction-0 line (35 of them) must be nonempty
   (this is E2'/Theorem T part 2 at w = 4 = d — machine-checked in M and
   re-checked here by the W kernel as its sanity anchor).
4. **Kill-restart demonstration (before the long run):** the runner writes a
   checkpoint every K = 400,000 supports (canonical JSON, sha256 over the
   canonical body, atomic tmp+rename+fsync); on restart it validates every
   ledger line (recompute sha256) and resumes at the exact next support
   index; a resume validator mode (`--validate-resume`) re-derives the
   count so far from the checkpoint and the frozen pre-statement's support
   ordering. DELIBERATE KILL: after the first checkpoint lands (≥400k
   supports processed), the operator kills -9 the runner, restarts it, and
   the resumed run must produce a count that CONVERGES to the same final
   total, with the ledger showing the restart boundary and zero
   re-counted supports (proven by per-slice counts and the final
   segmentation). This mirrors the delcap resume contract (validate_resume,
   checksum-validated lines, never truncate/repair; a checksum-invalid or
   gap-y ledger hard-aborts).

## 3. Support ordering (fixed, for the checkpoint ledger)

Supports of weight 4 over {0..139}, enumerated in LEXICOGRAPHIC order of
their sorted 4-tuples (itertools.combinations order). The checkpoint
records `next_index` = the number of supports fully processed (in this
order); a resumed run continues from exactly that index. The ordering is
frozen by THIS file (pre-compute); the runner asserts the count of the
first slice (index 0..399,999 resummarizes) matches a pre-computed offline
checksum of the 400,000th support tuple (recorded in the checkpoint's
first line: the tuple (0,1,2,3) at index 0, and a mid-index probe tuple
recorded with its index for the restart boundary audit).

## 4. Compute, resources, monitoring

- One single-core `nice -n 19` process, threads pinned to 1; the census is
  CPU-bound exact modular arithmetic (no GPU, no float).
- Pricing: C(140,4) = 15,329,615 supports; M's measured single-core rate at
  N=140 w≤3 was 603.5 supports/s (457,450 in 758.0 s); the w=4 per-support
  cost is higher (136-column rref per support vs the ≤3 window's average
  ~138.4 columns) so the honest estimate window is 2.7–5.5 h; the kill-
  restart demo adds ~5 min; per-support dims recorded for every nonempty
  support (expected count open — could be ~35 (lines only) or thousands;
  no committed prediction).
- Launch coordination: the launch WAITS for peer KgEvenBoxRerun's gate-B
  band to clear (its message: wall-clock-critical band [4.083,6.0], today);
  W runs nice -n 19 single-core, the kindest slot the machine gives; if the
  band does not clear before 2026-09-02T06:00Z, W's launch is postponed and
  the campaign hands back with the validated instrument + restart protocol
  (per assignment milestone (v)).
- Monitoring: the runner prints progress (index, count, rate) every 100k
  supports and on every checkpoint; the operator (this agent) polls the
  log; the checkpoint ledger carries per-slice times.
- Horizon: if the census will not finish within this session, the campaign
  freezes partial results + the resume protocol and hands back (Main
  re-launches a finishing agent). This is milestone (v); the kill-restart
  demo (milestone iii) is what makes a handback SAFE.

## 5. What this campaign does not claim (rule 7)

Swept: the weight-EXACTLY-4 point-support window at (421,(4,5,7),(1,1,1),
Λ=Id), complete census (15,329,615 supports), candidates = normalized
reduced basis of each nonzero V∩F^S, per-candidate delta EXACTLY closed
by exhaustive ascending-cost line-tuple sweep, ρ^window = min wt/delta
over nonzero supports (wt = 4 exactly). NOT swept: weights ≥ 5 (the w=4
result says NOTHING about whether heavier supports beat ρ^window_w4 —
ρ_inst stays open over the whole space); weights ≤ 3 (theorem-empty,
re-anchored completely by M's C-PN4 at the same instance — byte-copied
to the freeze, not re-run); every other (q,s,t,Λ); nonprime fields;
non-identity Λ; the global ρ_inst everywhere; Conjecture 4.2's content;
the TR26-150 paper claims (nothing here touches the paper's statements).
Λ-generality: inherited from M's L7 (pointwise-division bijection; the
emptiness/first-nonempty-weight statements are Λ-free; the RATIO ρ^window
is Λ-DEPENDENT in general and is only claimed at Λ=Id here).

## 6. Adjudication (pre-committed)

- Primary observed quantity: rho^window_w4 = min over the census's nonzero
  supports of Fraction(4, delta(S)) — EXACT (rationals; the delta sweep is
  exhaustive with the machine proof of absence of everything cheaper).
- Verdict vs H-GATE: H-GATE (as stated in P2, min-form) is ABOUT the ≤3
  window; the w=4 window is what the mechanism theorem predicts NONEMPTY
  (first nonempty weight = d = 4). The campaign therefore reports:
  (i) the theorem's expected signature — nonzero supports exist (at least
  the 35 lines) — with count and per-support dims; (ii) rho^window_w4
  exactly; (iii) the P-prediction status: no committed P prediction existed
  for this window (P's 3/5-iff-{2,3} was scoped to ≤3 windows; if the
  landed rho^window_w4 = 3/5, that is REMARKABLE-NEW evidence that P's
  3/5 signature reaches the w=4 window and would be reported as such, NOT
  as "P confirmed" — its scope never covered this window; if rho = 1,
  consistent with the ≤3-window history at non-{2,3} rows).
- Any support S with dim(V ∩ F^S) > 1 is recorded with its full reduced
  basis (multi-dim intersections are the "non-basis values" gap of gate-B's
  rule-7; exact delta closure is run per BASIS vector after normalization).
- Misses/failures of any prediction in §1 (P1-P3) are reported first and
  prominently, per the battery rule.

## 7. Option (b) — recorded follow-up candidate (per Main's instruction;
NOT run in this campaign)

A genuine min(s) ≤ 3 instance at N = 140 (e.g. q = 211, s = (2,7,10),
pairwise coprime... VERIFY coprime ⇒ s_i | q−1 = 210: 2|210, 7|210, 10|210
✓) — as a COMPARISON instance it would need its own paired censuses
(w ≤ 3 baseline + w = 4 window at the same instance) to attribute any
ratio movement to the window widening rather than the instance change;
priced as its own campaign with its own pre-statement. Recorded here as a
named candidate only: **"min3-N140 paired-window comparison"**.
