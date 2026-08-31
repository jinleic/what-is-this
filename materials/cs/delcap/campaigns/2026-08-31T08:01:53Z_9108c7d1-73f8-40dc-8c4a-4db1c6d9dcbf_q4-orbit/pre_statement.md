# delcap Pre-Statement — q=4 sparse orbit certificates and q=3 n=11 extension

Committed as the FIRST file in this campaign directory, before the q=4 population census, every startup anchor, and every new certificate row. Owner: DelcapNextQ. Date: 2026-08-31T08:01:53Z.

A pre-campaign development dry-run in `delcap/scratch/scope_bA.py` tested only the new exact orbit-sum identity on four tiny points `(q,n) in {(2,3),(3,3),(3,4),(2,5)}`. It produced no capacity row, no sandwich comparison, and no q=4 census. Its script and stdout will be frozen for provenance. The campaign computation starts only after this statement.

## Highest-value extension chosen

Candidate 1 from the ticket, the q=4 orbit ladder, is primary. Tavakoli-Nguyen-Bose report numerical experiments only for q in `{2,3}`; their Corollary 1 / Theorem 1 closed-form sandwich applies at q=4, so q=4 gives a new alphabet rather than two more rungs of an existing ladder. Candidate 2 is included only as the separately fixed secondary block `(q,n)=(3,11)` because it uses the identical instrument and was explicitly named mechanically reachable by the frozen q=3 campaign. Candidate 3, the cause of the TNB row-to-row rounding inconsistency, is not attacked: no author intermediates or code are available, so this instrument cannot ESTABLISH a cause (rule 14); the standing no-single-offset result is not reinterpreted or revived.

## Group, object and convention

The group is exactly `G = S_q x C_2` acting on SYMBOL VALUES, not positions: `sigma in S_q` relabels the alphabet and `rho` reverses the word; both act simultaneously on input and output. Words are tuples in left-to-right deletion order. The convention is

`A(x,y) = # { i_1 < ... < i_k : x[i_1]...x[i_k] = y }` and
`W(y|x) = A(x,y) d^(n-k) (1-d)^k`.

Deletion acts on POSITIONS. Therefore `A(sigma x, sigma y)=A(x,y)` and `A(rho x,rho y)=A(x,y)`. The falsified `S_n` input-POSITION/type reduction is excluded. Before any new row, an explicit position-permutation counterfactual must break in the same self-test that verifies the live value-permutation/reversal identities.

For an input orbit `O_j`, output orbit `T_t`, input representative `r_j`, and output representative `y_t`, the new implementation uses the exact identity

`F_j(y_t) = sum_{x in O_j} A(x,y_t)
            = |O_j|/|T_t| * sum_{y in T_t} A(r_j,y)`.

All terms are integers; divisibility is asserted. This avoids the frozen `merged_OS` sweep over all `q^n` input words while computing the SAME per-word orbit-constant marginal. It changes no certificate theorem. The identity is derived from orbit-stabilizer plus joint equivariance and must be equality-checked against the frozen per-word constructor before use.

## Stage 0 — exact q=4 population census, before any q=4 row

Exact counts are computed for q=4 at n=5..11:

* input words `4^n`;
* output words `sum_{k=0}^n 4^k`;
* input `S_4 x C_2` orbit representatives and flattened output-orbit representatives, by restricted-growth strings modulo reversal, with `sum orbit_size = 4^n` asserted;
* locator word support `|union_j supp A(r_j,.)|`;
* exact representative subsequence incidences `sum_j |supp A(r_j,.)|`;
* exact orbit-summed nonzero slots `sum_j |{T_t: sum_{y in T_t} A(r_j,y)>0}|`;
* dense-locator bytes `8 * n_input_orbits * word_support` and sparse-locator operation count per BA iteration.

Known-value anchors are asserted at startup, not afterwards: frozen `(q,n)=(3,10)` has exactly 5,002 input g-orbits and 7,566 flattened output g-orbits; q=4 restricted-growth counts equal `sum_{b=1}^4 S(n,b)`; orbit sizes sum to `4^n`; each channel row sums exactly to the common denominator. No q=4 certificate is run until the census artifact is written.

The census is the resource instrument, not a certificate verdict. Its fixed adjudication: n is called certifiable only if the population and one measured rehearsal fit 96 GiB RAM and 90 minutes per row under the frozen 4000-iteration locator. The row box below is NOT shrunk after observing capacity values. Any row that misses the resource stop is reported `NOT_REACHED_RESOURCE` with its counts and timing; it is never silently omitted.

## Startup gates before every new row

All must pass before the first q=4 or q=3,n=11 certificate:

1. Re-run the frozen orbit self-test A1-A6.
2. Load the frozen 20-row `orbit_rows.jsonl` and reproduce `(3,6,1/2)` and `(3,10,1/2)` with the frozen code path. `cert_lo_per_symbol`, `cert_hi_per_symbol`, `primal_ball`, `dual_ball`, and verdict must match bit-for-bit. Failure aborts.
3. Re-run the raw per-WORD snap guard at `(2,10,1/20)`, snap `2^-40`: exactly one reachable word mass is below resolution and snaps to zero; `cert_dual` must return `+infinity` (no claim), never a finite too-small upper bound.
4. Exhaustively verify q=4 symbol-value equivariance at `(q,n)=(4,3)`, including reversal. In the same path, the fixed position-permutation plant `x=0010 -> 0100`, `y=10` must give `A=1` versus `A=2`; if it does not break, the instrument cannot distinguish the falsified route and aborts.
5. Equality-check the exact orbit-sum construction against frozen `merged_OS` at `(3,9,1/2)` (all canonical output reps, all input reps; per-word samples; exact row-normalization closure). Run the new sparse locator/certificate path at `(3,6,1/2)` and `(3,10,1/2)` and require overlap with the frozen intervals.
6. Dense-vs-orbit cross-validation at q=4: `(4,5,d)` for every `d in {1/2,1/5,1/10,1/20}`, plus `(4,6,1/2)`. A disjoint pair falsifies the q=4 orbit implementation and aborts.
7. Permanent width-sign check on EVERY row: a finite dual upper endpoint must be at least the primal lower endpoint. Negative width is an invalid-certificate shape and aborts; it is never printed as a result.

## Fixed row box (rule 16)

Primary q=4 block, 24 rows:

* `(q,n) = (4,n)` for every `n in {5,6,7,8,9,10}`;
* every `d in {1/2,1/5,1/10,1/20}`.

Secondary q=3 block, 4 rows:

* `(q,n)=(3,11)`;
* every `d in {1/2,1/5,1/10,1/20}`.

Order is fixed: startup gates; q=4 ascending n and, within n, d in `(1/2,1/5,1/10,1/20)` order; then q=3,n=11 in the same d order. Resume may skip only rows already byte-present with a valid checksum. Nothing is added after a promising row. Rows at q=4,n=11 and q=3,n=12 are excluded even if a rehearsal looks cheap; they require a separate pre-statement.

## Certificates, precision and comparison

* Arb working precision: 400 bits, outward rounding.
* All probabilities entering the certificates: exact `fmpq` rationals.
* Input snap denominator `2^30`; output per-word snap denominator `2^30`; full-support `+1` bump and uniform fallback exactly as frozen.
* Locator: float64 BA, 4000 iterations, feeds only the rational snap. Its arithmetic is `COMPUTATIONAL-EVIDENCE` and has no bearing on certificate validity, only width.
* Primal: exact mutual information of the explicit rational input distribution; valid LOWER bound for any candidate.
* Dual: exact `max_x KL(W(.|x)||D')` over input g-orbit reps for an exact full-support G-invariant rational `D'`; equality of KL on each input orbit is an anchored identity. Any zero support on a reachable output returns `+infinity` (no claim).
* TNB comparison columns: exact-Arb recomputation of their Corollary 1 / Theorem 1 closed forms: `LB+ = (1-d)log2 q + H_Bin(n,1-d)/n - h2(d) + Delta_n(d)/n`, `UB=(1-d)log2 q`. No q=4 printed table value is invented.
* Each row records the certified interval, exact machine-reported width, certified margin `cert_lo - LB+` and `UB - cert_hi`, verdict, resource counts, and evidence label.

Adjudication, fixed before the first row:

* `CERT_LOWER_BEATS_LBplus` iff the Arb lower endpoint of `primal/n` is strictly above the Arb upper endpoint of `LB+`.
* `CERT_UPPER_BEATS_UB` iff the Arb upper endpoint of `dual/n` is strictly below the Arb lower endpoint of `UB`.
* Every miss, infinite dual, or resource stop is reported verbatim.
* If a certified interval lies outside the published sandwich rather than strictly inside it, stop and escalate to Main before any README statement. No tuning toward a published value.

## Budget and fixed resource stops

* Machine: Apple M3 Ultra, 96 GiB RAM, CPython 3.14.3, python-flint 0.9.0, NumPy 2.5.2.
* One sequential low-priority campaign process; no concurrent row jobs.
* Population census: hard stop 60 minutes; partial exact counts and the first uncounted rung are reported if exceeded.
* Per certificate row: hard stop 90 minutes. If any row reaches the stop, mark it `NOT_REACHED_RESOURCE` and continue only to smaller-or-equal pre-registered populations; never retarget the box.
* Whole certificate stage: 10 wall-clock hours; a budget exhaustion is an honest `FAILURE TO CERTIFY` for unfinished rows.
* No project-wide tests, formatters or linters.

## Rule-7 scope sentence, committed in advance

The intended sweep is exactly q=4 with `(n,d)` in `{5,6,7,8,9,10} x {1/2,1/5,1/10,1/20}`, plus q=3 with `(n,d)={11} x {1/2,1/5,1/10,1/20}`, using `G=S_q x C_2` on symbol VALUES and word reversal, Arb 400-bit outward rounding, exact rationals, input/output snaps `2^30`, and full finite channel support (no truncation, tail exactly zero). NOT swept: q>=5; q=4,n>=11; q=3,n>=12; other d; non-G-invariant input families as a separate search domain; the falsified S_n position/type reduction; TNB rounding cause; Morozov-Duman, Pinto-Ribeiro, and every n->infinity/asymptotic capacity claim. Every landed row is only a finite-n theorem.

## Anti-numerology and prospective/retrospective split

The population census, anchors, fixed grid, d order, arithmetic precision, snaps, verdict inequalities, and stop rules are prospective. Every row is reported including misses. Any hypothesis first noticed after rows land is labelled retrospective in the report and cannot change the box or adjudication. A counterfactual plant (position-permutation break and raw-snap +infinity branch) must produce its known nonzero/failure signal before a zero or pass is trusted.

Signed: DelcapNextQ. No q=4 census, startup anchor, or new capacity row has run as of this commit.