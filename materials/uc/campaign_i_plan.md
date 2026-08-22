# Campaign J plan — t = 0.3821660112501052 ≥ ψ + 2×10⁻⁴

**Status: LAUNCHED 2026-08-20 18:15 UTC** as campaign `certJ0–7`,
`uc/campaigns/cert3_20260820T181524Z_08ad738a1cf24a08bfc4189e65f9ef44_57a8908ba7b8/`,
code SHA-256 `57a8908ba7b8…` (the line-64 edit below makes it differ from
campaign I's `2f23a58ebdb8…`). This file's PREPARED sections are the record
of the plan as-executed: the launch conditions it demanded (t0 slice-7
COMPLETE on campaign H; no live workers contending) had both been superseded
comfortably by campaign I — which went 8/8 COMPLETE at t0 — by launch time.
The target is the rational `t = 0.3821660112501052`, which **strictly exceeds**
ψ + 2×10⁻⁴ (never equals it: t is rational, ψ is irrational) — see §2.

Everything in this file downstream of the exact-arithmetic check in §2 is
either machine-checked fact (cited to file+line or command output) or is
labeled **NUMERICAL** (sampled/projected evidence, not machine-checked).

---

## 1. Exact source edits needed (one line)

Baseline: `math/uc/cert3_par.py`, sha256
`a741aa167a79fb1058d89373db15e19d4124576a7eb7fb76c57ee3f70d3709bc`
(verified 2026-08-18; byte-identical to the copy frozen in campaign H's
snapshot, so the quoted lines are the same bytes H certified under). Fresh
read of the constants block `[math/uc/cert3_par.py#EE74]`:

```python
61: NSLICES = 8
62: MIN_WIDTH = 1e-3
63: WORK_PREC = 80
64: COLLAR_MIN_WIDTH = 1.25e-4
65: FACE_MIN_WIDTH = 1e-3
66: FACE_BOX_BUDGET = 100000
67: MAX_BOXES = 2_000_000_000
68: PROGRESS_EVERY = 4_000_000
```

### The edit (single line)

| line | current | change to | reason |
|---|---|---|---|
| 64 | `COLLAR_MIN_WIDTH = 1.25e-4` | `COLLAR_MIN_WIDTH = 6.25e-5` | one extra halving of the collar floor |

Mechanics: `COLLAR_MIN_WIDTH` flows from this constant into
`cert3.certify(..., collar_min_width=...)` at `cert3_par.py:575` and into
`launch.json` `parameters["collar_min_width"]` via `runtime_parameters`
(`cert3_par.py:241`). Per `cert3.py:1414-1416` (docstring): *"sink-collar-
touching boxes are refined to this finer width floor (sink faces carry the
global margin, but the corner bound's slop is first-order, so they need ~2-3
extra halvings)"*; the floor is applied per-coordinate at `cert3.py:1553-1554`
(`_floor(j)` returns `collar_min_width` for collar boxes, else `min_width`).

**NUMERICAL** rationale (NextTarget `runtime_risks_at_psi+2e-4`): the sampled
global min-Λ margin shrinks from ≈4.556×10⁻⁴ at ψ+10⁻⁴ to ≈2.926×10⁻⁴ at
ψ+2×10⁻⁴ (margin law 1.6298·(c\*−t), sampled evidence only). The old floor
1.25×10⁻⁴ would rise from 0.27× to 0.43× of the margin; the ~4.1 M
corner-cleared depth-42..47 terminals observed on G slice 7 hug that floor
and would need ~1 extra halving — if they strand as residuals the campaign
fails acceptance (`residual == 0` required). Halving to 6.25×10⁻⁵ = 0.21× of
the new margin restores ≈ the t0 ratio. Bounded extra cost: the ratio rule is
scale-free, so most collar boxes are unaffected.

### Lines deliberately NOT changed, with implications

| line | constant | value | implication at ψ+2×10⁻⁴ |
|---|---|---|---|
| 62 | `MIN_WIDTH` | `1e-3` | becomes 3.4× margin (was 2.2× at t0). Second-order centered rules need width ~ √margin, i.e. only ≈1.25× tighter than t0 (**NUMERICAL**). Keep. |
| 65 | `FACE_MIN_WIDTH` | `1e-3` | face margins track the same law; `face_bb` residual path (`cert3.py:1008-1010`) unchanged. Keep. |
| 66 | `FACE_BOX_BUDGET` | `100000` | no evidence of pressure (G slice 7: 20.2 M face nodes under this per-subtree budget). Keep. |
| 67 | `MAX_BOXES` | `2_000_000_000` | slice-7 honest bracket 44 M–530 M boxes, worst joint corner 530 M < 2×10⁹ (**NUMERICAL**); binds only at ψ+3×10⁻⁴ pessimistic scenarios. Keep. |

Do **not** touch `cert3.py` (its internal default `collar_min_width = 1.25e-4`
at `cert3.py:1440-1442` is dead for campaign runs — `cert3_par.py:575` always
passes the constant explicitly). Editing any `EXECUTABLE_FILE`
(`cert3_par.py` at line 64 included) changes `code_sha256`
(`combined_source_hash`), which is by-design: campaign I gets a fresh id and
a fresh immutable snapshot at `create` time; H's frozen snapshot and verdict
are unaffected.

### Time budget: 259200 s

The time budget is **not** a source constant — it is the `SECONDS` argument to
`create` (`cert3_par.py:446-449, 730-735`), recorded per run as
`time_budget_s` (`cert3_par.py:477`) and validated as a positive integer
(`cert3_par.py:393-395`). Campaign I uses **259200** (= 72 h), satisfying the
NextTarget gate "per-slice time_budget must be ≥ 72 h for a 3-day target
(COMPLETE requires `budget_time == 0`)". (H used 86400 s = 24 h.)

---

## 2. Exact target check — VERIFIED (machine-checked, 2026-08-18)

`certified_t_decimal(offset)` (`cert3_par.py:213-224`) walks floats upward
from `float(ψ+offset)` and returns the first `repr` for which
`target_relation_holds` (`cert3_par.py:198-210`) holds as **exact rational
arithmetic**: with `r = Fraction(t_decimal) − Fraction(offset)`, require
`0 ≤ r ≤ 3/2` and `r² − 3r + 1 ≤ 0` (r at least the smaller root ψ of
x²−3x+1, i.e. t ≥ ψ+offset exactly).

Note: the certificate's `t` is a **rational** such that `t ≥ ψ + 2×10⁻⁴`;
it never equals ψ+2×10⁻⁴ pointwise, since ψ is irrational. The printed
`t_decimal 0.3821660112501052` (below) satisfies `r = t − 0.0002 =
0.3819660112501052 > ψ` and `r² − 3r + 1 < 0` strictly — hence
`t > ψ + 2×10⁻⁴` in exact rational arithmetic.

Command run (the single permitted check, `.venv` python, < 0.1 s):

```sh
cd /Users/jinleic/jinleic-workspace/math && ./.venv/bin/python -c "
from fractions import Fraction
from decimal import Decimal
import math
def target_relation_holds(t_decimal, offset):
    t = Fraction(Decimal(t_decimal)); o = Fraction(Decimal(offset))
    r = t - o
    if r < 0 or r > Fraction(3,2): return False
    return r*r - 3*r + 1 <= 0
print('holds:', target_relation_holds('0.3821660112501052', '0.0002'))
r = Fraction(Decimal('0.3821660112501052')) - Fraction(Decimal('0.0002'))
print('r =', r, '=', float(r))
print('poly r^2-3r+1 =', r*r - 3*r + 1)
"
```

Observed output:

```
holds: True
r = 954915028125263/2500000000000000 = 0.3819660112501052
poly r^2-3r+1 = -673679581180831/6250000000000000000000000000000
```

Interpretation: r = 0.3821660112501052 − 0.0002 = 0.3819660112501052 exactly,
and r²−3r+1 < 0 **strictly**, so r > ψ = (3−√5)/2 exactly. Therefore
t = 0.3821660112501052 ≥ ψ + 2×10⁻⁴ holds in exact rational arithmetic.

Expected consequence: `create 0.0002 259200` records
`t_decimal == "0.3821660112501052"` (the next-lower float's decimal
0.38216601125010515 < ψ+2×10⁻⁴ fails the relation, so the upward walk lands
here). **Defensive gate:** if the `create` output prints any other
`t_decimal`, do not launch workers — rerun the one-liner above on the printed
value and only proceed if it prints `holds: True`. Every worker and the
collector re-derive and re-check this internally anyway
(`verify_launch` → `runtime_parameters(offset)`, `cert3_par.py:326-327`,
`432-439`), so an inadmissible t cannot certify.

---

## 3. Exact launch sequence

**As executed 2026-08-20 18:15 UTC.** The original gating (wait for campaign
H's frozen collector to print `COMPOSITE CERTIFICATE` before creating the
ψ + 2×10⁻⁴ campaign) was superseded: by the time of launch the stronger
campaign I had already committed all 8 slices at t0 with zero residual, and
the only live processes were I's frozen collector's own replay on ~1 core —
the machine had ~27 free cores.

**Once — create the campaign (this is the only step that writes a new
campaign dir; it also applies nothing to any existing campaign):**

```sh
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -B uc/cert3_par.py create 0.0002 259200
```

This was executed 2026-08-20 18:15 UTC and printed campaign path
`uc/campaigns/cert3_20260820T181524Z_08ad738a1cf24a08bfc4189e65f9ef44_57a8908ba7b8/`,
`CODE_SHA256 57a8908ba7b8…` (≠ H/I's `2f23a58ebdb8…` because of the line-64
edit) and eight `RUN` lines, then workers `certJ0–7` were launched from the
frozen snapshot's `cert3_par.py`, logs outside the campaign dir.

**Once — capture the fresh paths (substitute the LAUNCH_PATH just printed):**

```sh
LP=/Users/jinleic/jinleic-workspace/math/uc/campaigns/cert3_<NEW_ID>/launch.json
DRV=${LP%/launch.json}/snapshot/cert3_par.py
```

**Then two paste commands per slice** (slice i = 0…7; all eight may run
concurrently; each supervisor spawns one single-thread worker):

```sh
cd /Users/jinleic/jinleic-workspace/math
nohup ./.venv/bin/python -B "$DRV" run "$LP" 0 > /tmp/certI_slice0.log 2>&1 &
```
```sh
cd /Users/jinleic/jinleic-workspace/math
nohup ./.venv/bin/python -B "$DRV" run "$LP" 1 > /tmp/certI_slice1.log 2>&1 &
```
```sh
cd /Users/jinleic/jinleic-workspace/math
nohup ./.venv/bin/python -B "$DRV" run "$LP" 2 > /tmp/certI_slice2.log 2>&1 &
```
```sh
cd /Users/jinleic/jinleic-workspace/math
nohup ./.venv/bin/python -B "$DRV" run "$LP" 3 > /tmp/certI_slice3.log 2>&1 &
```
```sh
cd /Users/jinleic/jinleic-workspace/math
nohup ./.venv/bin/python -B "$DRV" run "$LP" 4 > /tmp/certI_slice4.log 2>&1 &
```
```sh
cd /Users/jinleic/jinleic-workspace/math
nohup ./.venv/bin/python -B "$DRV" run "$LP" 5 > /tmp/certI_slice5.log 2>&1 &
```
```sh
cd /Users/jinleic/jinleic-workspace/math
nohup ./.venv/bin/python -B "$DRV" run "$LP" 6 > /tmp/certI_slice6.log 2>&1 &
```
```sh
cd /Users/jinleic/jinleic-workspace/math
nohup ./.venv/bin/python -B "$DRV" run "$LP" 7 > /tmp/certI_slice7.log 2>&1 &
```

Protocol invariants baked into the tooling (no operator discipline needed):
`-B` everywhere so no `__pycache__` pollutes the pristine snapshot (the
listing gate at `cert3_par.py:410-420` rejects any extra entry forever);
logs go to `/tmp`, not the campaign dir; supervisors refuse to overwrite
committed results/traces (`cert3_par.py:670-675`); the worker re-verifies the
frozen bytes before committing (`cert3_par.py:588`). Never `import` snapshot
modules from another program — run them only as scripts.

**After all eight slices commit** (each prints `RESULT_PATH`), run the frozen
collector (same line `create` printed as `COLLECT`):

```sh
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -B "${LP%/launch.json}/snapshot/cert3_collect.py" "$LP"
```

Exit 0 with `COMPOSITE CERTIFICATE` is the only acceptance (schema v2
collector replays all eight traces itself). Post-launch bookkeeping (not part
of this plan's paste block): add the campaign row to
`uc/campaigns/README.md`, and an `run_i_watch.py`-style watcher plus
`PROOF.md` target bump only after certification.

---

## 4. Rollback note

**Nothing to roll back — the campaign is live.** The campaign J directory
was created fresh by §3's `create` on 2026-08-20 18:15 UTC; every earlier
campaign's directory and snapshot are immutable and untouched by this plan.
The only working-tree change is the one-line constant at `cert3_par.py:64`;
since `create` froze the new snapshot, later working-tree edits cannot affect
the running campaign (workers verify against their own snapshot copy,
`cert3_par.py:408-409`). If campaign J is abandoned mid-flight: leave its
directory in place as immutable evidence, record the verdict
(`NO CERTIFICATE` / `SUPERSEDED`) in `uc/campaigns/README.md`, and revert
line 64 to `1.25e-4` only if a future campaign wants the old floor. No
earlier-campaign artifact is ever at risk.

---

## 5. Expected runtime bands — all NUMERICAL (NextTarget analysis)

Projections from the margin law 1.6298·(c\*−t) (sampled minima reproduce to
4 digits), cost exponent α∈[2,4] (structural argument from G's rule mix; no
measured t-sweep), and H's measured splitter gains (278.7× on the degenerate
residual slab, 2.0–2.8× slab throughput) — none machine-checked.

| quantity | NUMERICAL projection |
|---|---|
| slice-7 tree at ψ+2×10⁻⁴ | central **~110 M boxes** (0.08–0.15 B); honest bracket 44 M–530 M |
| slice-7 wall time | **25–38 h loaded** (800–1,200 boxes/s); worst-joint corner 54–184 h |
| time budget | 259200 s = 72 h per slice — covers the central case with ≈2× headroom; worst-joint corner (184 h) would exceed it and end `INCOMPLETE` with `budget_time > 0` |
| slices 0–6 | start 2–3.5× smaller than slice 7 at t0 and scale by the same multiplier ⇒ makespan = slice 7 alone |
| aggregate | ~330–700 M box-events total vs ~1,140 M capacity (8 workers × 3 days × ~55 % core share × ~1,000 b/s) |
| live rate anchor | H slice 7: ~1,075 b/s loaded (18.74 M trace bytes / 4.84 h, single sample) |
| hard ceiling context | ψ+3×10⁻⁴ needs ~0.9–2.1 B boxes on slice 7 (215–323 h) — does not fit 72 h; ψ+4×10⁻⁴ exceeds c\* ≈ ψ+3.795×10⁻⁴, mathematically void |

Failure modes to watch (NUMERICAL): corner-debris stranding at the collar
floor (mitigated by the §1 edit — the point of it), and external load
(load-avg ~100 on 28 cores; a dedicated core per worker would ~2× all wall
times). Acceptance still requires the machine-checked tallies `residual == 0`,
`stack == 0`, `budget_boxes == 0`, `budget_time == 0` on all eight slices
plus eight independent trace replays inside the collector.
