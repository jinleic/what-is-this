#!/usr/bin/env bash
#
# closeout.sh — one-shot closeout for campaign I
# (cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8).
#
# Run ONCE, after campaign I's frozen collector has printed
# COMPOSITE CERTIFICATE into the campaign directory's collect_output.txt
# (the iwatch watcher writes it when all eight slices replay PASS).
#
# What it does, in order (nothing is written before every check has passed):
#   1. Refuses to run unless collect_output.txt contains the ANCHORED success
#      marker '^COMPOSITE CERTIFICATE:' — a substring match would also accept
#      the collector's failure line 'NO COMPOSITE CERTIFICATE -- N problem(s):'.
#   2. Cheap structural validation of that SAME output (no re-execution):
#      the collector is the verification gate — it already replayed all 8
#      traces in-process before printing.  Third-party evidence would use a
#      fresh collector run, not this closeout script; replaying 488 M nodes a
#      second time here proves nothing new and costs hours.  This step checks
#      the output's structure: CHECKER_SHA256 (64-hex), 'ALL MACHINE CHECKS
#      PASS', '^COMPOSITE CERTIFICATE:', exactly 8 'REPLAY slice=i PASS'
#      lines, and per-slice processed tallies equal to the committed records.
#   3. Preflight of EVERY ledger target state, read-only — PROGRESS.md's
#      insertion anchor exists, the campaigns/README.md campaign-I row is
#      exactly one row in its pre-CERTIFIED pending form, uc/README.md has a
#      level-1 heading, and CERTIFICATE.txt is either absent or already
#      byte-identical to the collector output.  Any deviation → die before
#      touching anything.
#   4. Only then, archive: dumps the collector output verbatim to
#      uc/campaigns/CERTIFICATE.txt (NO-REPLACE: identical bytes → no-op,
#      different pre-existing bytes → refuses loudly).
#   5. Only after validation + preflight + archive pass does it mutate any
#      record: appends the certificate paragraph to math/PROGRESS.md under a
#      new dated session-9 heading, flips the campaign-I pending row in
#      math/uc/campaigns/README.md to CERTIFIED (date + slice table), and
#      updates the headline of math/uc/README.md (certificate exists at
#      t = psi + 1e-4, campaign id, collector SHA-256).
#
# Idempotent: every mutation is guarded by the sentinel comment
#   <!-- closeout:cert3_20260818T212601Z_425f109c -->
# so rerunning appends nothing twice and the structural check is cheap.
#
# No campaign artifact, snapshot file, or cert*.py source is ever written.

set -euo pipefail
umask 022

UC_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MATH_DIR="$(cd -- "$UC_DIR/.." && pwd)"

CAMPAIGN_ID="20260818T212601Z_425f109c15b64a6198785c6cebbbdaab"
CODE_HASH="2f23a58ebdb8"
CAMPAIGN_DIR="$UC_DIR/campaigns/cert3_${CAMPAIGN_ID}_${CODE_HASH}"
COLLECT_OUT="$CAMPAIGN_DIR/collect_output.txt"
COLLECTOR="$CAMPAIGN_DIR/snapshot/cert3_collect.py"
LAUNCH_JSON="$CAMPAIGN_DIR/launch.json"
CERT_TXT="$UC_DIR/campaigns/CERTIFICATE.txt"

SENTINEL='<!-- closeout:cert3_20260818T212601Z_425f109c -->'
PY="$MATH_DIR/.venv/bin/python"

DATE="$(date -u +%F)"
DATE_FULL="$(date -u '+%Y-%m-%d %H:%M UTC')"

say() { printf '[closeout %s] %s\n' "$(date -u +%H:%M:%S)" "$*"; }
die() { printf '[closeout] FATAL: %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Step 1 — gate: refuse unless the frozen collector already certified.
# ---------------------------------------------------------------------------
[ -x "$PY" ] || die "venv python not executable: $PY"
[ -f "$COLLECTOR" ] || die "frozen collector missing: $COLLECTOR"
[ -f "$LAUNCH_JSON" ] || die "launch manifest missing: $LAUNCH_JSON"

if [ ! -f "$COLLECT_OUT" ] || ! grep -q '^COMPOSITE CERTIFICATE:' "$COLLECT_OUT"; then
    die "refusing to run: $COLLECT_OUT lacks the anchored success marker '^COMPOSITE CERTIFICATE:'. (A substring match would also accept the failure line 'NO COMPOSITE CERTIFICATE -- ...'.) Campaign I is not certified yet."
fi
say "gate passed: anchored success marker '^COMPOSITE CERTIFICATE:' present (campaign $CAMPAIGN_ID)"

# ---------------------------------------------------------------------------
# Step 2 — structural validation of the SAME collector output (no replay).
#          The collector is the verification gate: it replayed all eight
#          traces in-process before printing.  A second 488 M-node replay
#          duplicates hours of work and proves nothing beyond the exit-0,
#          CHECKER-pinned output already on disk — so closeout.sh validates
#          that output's structure against the committed records instead.
# ---------------------------------------------------------------------------
export CLOSEOUT_MATH_DIR="$MATH_DIR" CLOSEOUT_UC_DIR="$UC_DIR" \
       CLOSEOUT_CAMPAIGN_DIR="$CAMPAIGN_DIR" CLOSEOUT_DATE="$DATE" \
       CLOSEOUT_DATE_FULL="$DATE_FULL" CLOSEOUT_SENTINEL="$SENTINEL"
if ! "$PY" - "$COLLECT_OUT" <<'PY_SUMMARY'
import json, os, re, sys

MATH = os.environ["CLOSEOUT_MATH_DIR"]
CAMP = os.environ["CLOSEOUT_CAMPAIGN_DIR"]
DATE_FULL = os.environ["CLOSEOUT_DATE_FULL"]
NSLICES = 8
H_REL = "cert3_20260817T212001Z_4e6268eb7ce44aa3883ef15ea401c7c7_2f23a58ebdb8"

def die(msg):
    print("[closeout] FATAL: %s" % msg, file=sys.stderr)
    sys.exit(1)

text = open(sys.argv[1], encoding="utf-8").read()
lines = text.splitlines()

# --- validate the archived collector output ------------------------------
checker = None
for line in lines:
    if line.startswith("CHECKER_SHA256 "):
        checker = line.split(None, 1)[1].strip()
if not checker or not re.fullmatch(r"[0-9a-f]{64}", checker):
    die("collector output lacks a valid CHECKER_SHA256 line")
if "ALL MACHINE CHECKS PASS:" not in lines:
    die("collector output lacks 'ALL MACHINE CHECKS PASS:'")
if not any(l.startswith("COMPOSITE CERTIFICATE") for l in lines):
    die("collector output lacks 'COMPOSITE CERTIFICATE'")

replayed = {}
for line in lines:
    m = re.match(r"^REPLAY slice=(\d+) PASS processed=(\d+) split=(\d+) "
                 r"face=(\d+) infeasible=(\d+) elapsed=([0-9.]+)s$", line)
    if m:
        replayed[int(m.group(1))] = {"processed": int(m.group(2)),
                                     "elapsed_s": float(m.group(6))}
if sorted(replayed) != list(range(NSLICES)):
    die("collector output does not carry exactly 8 REPLAY PASS lines")

# --- committed records --------------------------------------------------
launch = json.load(open(os.path.join(CAMP, "launch.json"), encoding="utf-8"))
runs = {r["slice"]: r for r in launch["runs"]}
if sorted(runs) != list(range(NSLICES)):
    die("launch.json does not name exactly slices 0..7")
rows = []
total = 0
for i in range(NSLICES):
    rec = json.load(open(os.path.join(CAMP, runs[i]["result_filename"]),
                         encoding="utf-8"))
    t = rec["tallies"]
    if rec.get("verdict") != "COMPLETE" or rec.get("worker_exit_code") != 0:
        die("slice %d record is not COMPLETE with worker exit 0" % i)
    if t["processed"] != replayed[i]["processed"]:
        die("slice %d: collector replay processed %d != committed %d"
            % (i, replayed[i]["processed"], t["processed"]))
    rows.append({"slice": i,
                 "w": "[%s,%s]" % (runs[i]["w_lo"], runs[i]["w_hi"]),
                 "processed": t["processed"],
                 "residual": t["residual"],
                 "stack": t["stack"],
                 "elapsed_s": replayed[i]["elapsed_s"]})
    total += t["processed"]
replay_total = sum(r["elapsed_s"] for r in rows)

# --- determinism witness vs campaign H ----------------------------------
h_dir = os.path.join(os.path.dirname(CAMP), H_REL)
same = []
if os.path.isdir(h_dir):
    hlaunch = json.load(open(os.path.join(h_dir, "launch.json"), encoding="utf-8"))
    hruns = {r["slice"]: r for r in hlaunch["runs"]}
    for i in range(NSLICES):
        hp = os.path.join(h_dir, hruns[i]["result_filename"])
        if os.path.exists(hp):
            hrec = json.load(open(hp, encoding="utf-8"))
            if hrec.get("trace_sha256") == json.load(
                    open(os.path.join(CAMP, runs[i]["result_filename"]),
                         encoding="utf-8")).get("trace_sha256"):
                same.append(i)

bar = "=" * 78
print(bar)
print("CAMPAIGN I CLOSEOUT — FINAL SELF-CHECK SUMMARY  (%s)" % DATE_FULL)
print("collector : %s" % os.path.join("snapshot", "cert3_collect.py"))
print("checker   : SHA-256 %s" % checker)
print("campaign  : %s" % launch["campaign_id"])
print("exit code : 0 (iwatch)  [ALL MACHINE CHECKS PASS] [COMPOSITE CERTIFICATE]")
print("-" * 78)
print("slice  w-interval     processed    residual  stack   replay elapsed")
for r in rows:
    print("  %d    %-13s %11s  %8s  %5s  %11s"
          % (r["slice"], r["w"], format(r["processed"], ","),
             format(r["residual"], ","), r["stack"],
             format(int(round(r["elapsed_s"])), ",") + " s"))
print("-" * 78)
print("total boxes      : %s   (8/8 slices COMPLETE; replayed tallies exact)"
      % format(total, ","))
print("replay expense   : %.1f h total inside the collector (as printed by it)"
      % (replay_total / 3600.0))
if same:
    print("determinism      : slices %s trace SHA-256 identical to campaign H"
          % ",".join(map(str, same)))
elif os.path.isdir(h_dir):
    print("determinism      : NO slice trace matches campaign H (unexpected)")
print(bar)
PY_SUMMARY
then
    die "collector output failed validation; no ledger file was touched"
fi

# ---------------------------------------------------------------------------
# Step 3 — preflight every ledger target state (read-only).
# ---------------------------------------------------------------------------
say "validation passed: collector output structure and tallies match the committed records"

"$PY" - <<'PY_PREFLIGHT'
import json, os, re, sys

MATH = os.environ["CLOSEOUT_MATH_DIR"]
UC = os.environ["CLOSEOUT_UC_DIR"]
CAMP = os.environ["CLOSEOUT_CAMPAIGN_DIR"]
SENT = os.environ["CLOSEOUT_SENTINEL"]
PROG = os.path.join(MATH, "PROGRESS.md")
CREADME = os.path.join(UC, "campaigns", "README.md")
UREADME = os.path.join(UC, "README.md")
CERT_TXT = os.path.join(UC, "campaigns", "CERTIFICATE.txt")
COLLECT_OUT = os.path.join(CAMP, "collect_output.txt")
ROW_PREFIX = "| 2026-08-18 21:26 |"

def die(msg):
    print("[closeout] FATAL: %s" % msg, file=sys.stderr)
    sys.exit(1)

# PROGRESS.md: if the sentinel is absent, the insertion anchor must exist.
cur = open(PROG, encoding="utf-8").read()
if SENT not in cur and "\n## " not in cur:
    die("PROGRESS.md has no '## ' heading to insert before")

# campaigns/README.md: if the sentinel is absent, exactly one campaign-I row
# in a recognised pre-CERTIFIED form must be present.
cur = open(CREADME, encoding="utf-8").read()
if SENT not in cur:
    rows = [l for l in cur.splitlines() if l.startswith(ROW_PREFIX)]
    if len(rows) != 1:
        die("expected exactly one campaign-I row in campaigns/README.md, found %d" % len(rows))
    if ("LIVE" not in rows[0]) and ("ALL 8 SLICES COMPLETE" not in rows[0]):
        die("campaign-I row is in an unrecognised state: %r" % rows[0][:80])

# uc/README.md: if the sentinel is absent, a level-1 heading must exist.
lines = open(UREADME, encoding="utf-8").read().splitlines()
if SENT not in "\n".join(lines) and not any(l.startswith("# ") for l in lines):
    die("uc/README.md has no level-1 heading")

# CERTIFICATE.txt: absent, or byte-identical to the collector output.
if os.path.exists(CERT_TXT):
    if open(CERT_TXT, "rb").read() != open(COLLECT_OUT, "rb").read():
        die("CERTIFICATE.txt exists with DIFFERENT bytes — refusing to overwrite (no-replace policy)")

print("[closeout] preflight passed: all ledger targets are in their expected pre-states")
PY_PREFLIGHT
[ $? -eq 0 ] || die "preflight failed — no ledger file was touched"

# ---------------------------------------------------------------------------
# Step 4 — dump the collector output to uc/campaigns/CERTIFICATE.txt.
#          NO-REPLACE: if the archive already exists it must be byte-identical
#          to the collector output (idempotent rerun → silent no-op); any
#          other case (different pre-existing certificate) refuses loudly —
#          though preflight already excluded it.
# ---------------------------------------------------------------------------
if [ -e "$CERT_TXT" ]; then
    if cmp -s "$COLLECT_OUT" "$CERT_TXT"; then
        say "CERTIFICATE.txt already present and byte-identical — archive step is a no-op"
    else
        die "CERTIFICATE.txt exists with DIFFERENT bytes — refusing to overwrite (no-replace policy)"
    fi
else
    cp "$COLLECT_OUT" "$CERT_TXT"
    say "collector output archived to uc/campaigns/CERTIFICATE.txt ($(wc -c < "$CERT_TXT" | tr -d ' ') bytes)"
fi

# ---------------------------------------------------------------------------
# Step 5 — ledger + README updates (guarded by the sentinel; idempotent).
# ---------------------------------------------------------------------------
"$PY" - <<'PY_MUTATE'
import json, os, sys

MATH = os.environ["CLOSEOUT_MATH_DIR"]
UC = os.environ["CLOSEOUT_UC_DIR"]
CAMP = os.environ["CLOSEOUT_CAMPAIGN_DIR"]
DATE = os.environ["CLOSEOUT_DATE"]
DATE_FULL = os.environ["CLOSEOUT_DATE_FULL"]
SENT = os.environ["CLOSEOUT_SENTINEL"]

NSLICES = 8
H_REL = "cert3_20260817T212001Z_4e6268eb7ce44aa3883ef15ea401c7c7_2f23a58ebdb8"
PROG = os.path.join(MATH, "PROGRESS.md")
CREADME = os.path.join(UC, "campaigns", "README.md")
UREADME = os.path.join(UC, "README.md")
RULES = ("infeasible", "corner", "center", "center_mixed", "center_mixed_swap",
         "center_w", "ratio", "face", "split")

def die(msg):
    print("[closeout] FATAL: %s" % msg, file=sys.stderr)
    sys.exit(1)

def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)

def rel(path):
    return os.path.relpath(path, MATH)

def write(path, text):
    tmp = path + ".closeout.tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(tmp, path)

launch = load(os.path.join(CAMP, "launch.json"))
runs = {r["slice"]: r for r in launch["runs"]}
if sorted(runs) != list(range(NSLICES)):
    die("launch.json does not name exactly slices 0..7")

rows = []
total = 0
for i in range(NSLICES):
    rec = load(os.path.join(CAMP, runs[i]["result_filename"]))
    t = rec["tallies"]
    if rec.get("verdict") != "COMPLETE" or rec.get("worker_exit_code") != 0:
        die("slice %d record is not COMPLETE with worker exit 0" % i)
    for k in ("stack", "residual", "budget_boxes", "budget_time"):
        if t.get(k) != 0:
            die("slice %d tally %s = %r, expected 0" % (i, k, t.get(k)))
    if t["processed"] <= 0 or rec.get("trace_size") != t["processed"]:
        die("slice %d processed/trace_size mismatch" % i)
    rows.append({"slice": i,
                 "w": "[%s,%s]" % (runs[i]["w_lo"], runs[i]["w_hi"]),
                 "t": t,
                 "sha": rec["trace_sha256"],
                 "elapsed_s": int(round(t["elapsed_ms"] / 1000.0))})
    total += t["processed"]

# --- archived collector output ------------------------------------------
co_lines = open(os.path.join(CAMP, "collect_output.txt"),
                encoding="utf-8").read().splitlines()
checker = None
cert = []
for line in co_lines:
    if line.startswith("CHECKER_SHA256 "):
        checker = line.split(None, 1)[1].strip()
    if line.startswith(("COMPOSITE CERTIFICATE", "at certified rational",
                        "The proved orbit-swap")):
        cert.append(line)
if not checker or len(checker) != 64 or any(c not in "0123456789abcdef" for c in checker):
    die("collect_output.txt lacks a CHECKER_SHA256 line")
if len(cert) < 3:
    die("collect_output.txt lacks the certificate sentence lines")
cert = cert[:3]
if launch["t_decimal"] not in cert[1] or launch["level_offset"] not in cert[1]:
    die("certificate sentence disagrees with launch.json t/offset")

# --- determinism witness vs campaign H -----------------------------------
h_dir = os.path.join(UC, "campaigns", H_REL)
same, hdiff_complete = [], []
if os.path.isdir(h_dir):
    hruns = {r["slice"]: r for r in load(os.path.join(h_dir, "launch.json"))["runs"]}
    for i in range(NSLICES):
        hp = os.path.join(h_dir, hruns[i]["result_filename"])
        if not os.path.exists(hp):
            continue
        hrec = load(hp)
        if hrec.get("trace_sha256") == rows[i]["sha"]:
            same.append(i)
        elif hrec.get("verdict") == "COMPLETE":
            hdiff_complete.append(i)

fmt = lambda n: format(n, ",")

# ========================================================================
# 4a. math/PROGRESS.md — new dated session-9 heading + certificate paragraph
# ========================================================================
L = []
L.append("## %s — Session 9 closeout: CAMPAIGN I CERTIFIED — COMPOSITE CERTIFICATE at certified rational t = %s ≥ ψ + 10⁻⁴" % (DATE, launch["t_decimal"]))
L.append(SENT)
L.append("")
L.append("Campaign `uc/campaigns/cert3_%s_%s/` (schema v2, 72 h budget per slice,"
         % (launch["campaign_id"], launch["code_sha256"][:12]))
L.append("executable code hash `%s…%s`): the frozen collector" %
         (launch["code_sha256"][:12], launch["code_sha256"][-6:]))
L.append("`snapshot/cert3_collect.py` re-proved all eight committed traces and printed")
L.append("**COMPOSITE CERTIFICATE**. Checker SHA-256 `%s`; launch SHA-256" % checker)
L.append("`%s`. Verbatim collector output archived at" % launch["launch_sha256"])
L.append("[`uc/campaigns/CERTIFICATE.txt`](uc/campaigns/CERTIFICATE.txt) (copy of the")
L.append("campaign's `collect_output.txt`, written by the `iwatch` replay watcher).")
L.append("")
L.append("The certified statement, verbatim from the collector:")
L.append("")
for c in cert:
    L.append("> %s" % c)
L.append("")
L.append("All eight slices COMPLETE with `worker_exit_code` 0 and")
L.append("stack = residual = budget_boxes = budget_time = 0; every trace's SHA-256,")
L.append("size-equals-processed, and full mathematical replay validated by the")
L.append("collector. Slice tallies as committed in the eight result JSONs:")
L.append("")
L.append("| slice | w-interval | processed | residual | stack | elapsed | trace SHA-256 |")
L.append("|---|---|---|---|---|---|---|")
for r in rows:
    L.append("| %d | %s | %s | 0 | 0 | %s s | `%s` |"
             % (r["slice"], r["w"], fmt(r["t"]["processed"]),
                fmt(r["elapsed_s"]), r["sha"][:16]))
L.append("")
L.append("Discharge-rule mix (same records):")
L.append("")
L.append("| slice | infeasible | corner | center | mixed | swap | w-only | ratio | face | split |")
L.append("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    t = r["t"]
    L.append("| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |"
             % (r["slice"], fmt(t["infeasible"]), fmt(t["corner"]),
                fmt(t["center"]), fmt(t["center_mixed"]),
                fmt(t["center_mixed_swap"]), fmt(t["center_w"]),
                fmt(t["ratio"]), fmt(t["face"]), fmt(t["split"])))
L.append("")
L.append("**Total %s boxes**, zero residual leaves anywhere." % fmt(total))
L.append("")
if same:
    L.append("Determinism witness: trace SHA-256 identical to campaign H's")
    L.append("(`cert3_20260817T212001Z_…_%s`, same code hash) on slices %s; the" %
             (launch["code_sha256"][:12], ", ".join(map(str, same))))
    L.append("remaining slices have no COMPLETE H counterpart (H was wall-clock")
    L.append("limited there), so no SHA equality is expected on them.")
    if hdiff_complete:
        L.append("**UNEXPECTED**: slices %s differ from H's COMPLETE records —"
                 % ", ".join(map(str, hdiff_complete)))
        L.append("investigate before announcing.")
elif os.path.isdir(h_dir):
    L.append("Determinism note: no slice trace SHA-256 matches campaign H's")
    L.append("(unexpected for the same code hash — investigate before announcing).")
L.append("")
L.append("Closeout validation (`uc/closeout.sh`, %s): the archived collector" % DATE_FULL)
L.append("output was structurally re-checked BEFORE this entry was written — the")
L.append("`COMPOSITE CERTIFICATE` marker, its CHECKER hash print, and all eight")
L.append("`REPLAY slice=i PASS` lines, with processed tallies equal to the")
L.append("records above.  No second full replay was executed: the collector's own")
L.append("in-process replay of all eight traces is the verification gate; a")
L.append("cheap structural re-check adds the bookkeeping guarantee without")
L.append("duplicating 488 M nodes of work.")
L.append("")
section = "\n".join(L)

cur = open(PROG, encoding="utf-8").read()
if SENT in cur:
    print("[closeout] %s already carries the campaign-I entry — skipped" % rel(PROG))
else:
    idx = cur.find("\n## ")
    if idx < 0:
        die("PROGRESS.md has no '## ' heading to insert before")
    new = cur[:idx + 1] + section + "\n\n" + cur[idx + 1:]
    write(PROG, new)
    print("[closeout] %s — session-9 closeout entry inserted" % rel(PROG))

# ========================================================================
# 4b. uc/campaigns/README.md — flip the LIVE row for campaign I to CERTIFIED
# ========================================================================
row_prefix = "| 2026-08-18 21:26 |"
inline = " / ".join(fmt(r["t"]["processed"]) for r in rows)
new_row = ("| 2026-08-18 21:26 | `…%s` | %s | v2 | 72 h | 8 results + 8 traces | "
           "**CERTIFIED %s** — frozen collector printed `COMPOSITE CERTIFICATE` "
           "(all 8 traces replayed in its own gate; closeout.sh structurally "
           "re-validated that output against the committed records): slices "
           "0–7 = %s boxes, **%s total**, "
           "residual 0 · stack 0 · `budget_time` 0 on every slice; checker "
           "SHA-256 `%s…`; verbatim output at [`CERTIFICATE.txt`](CERTIFICATE.txt) %s |"
           % (launch["code_sha256"][:12], launch["code_sha256"][:12], DATE,
              inline, fmt(total), checker[:16], SENT))

cur = open(CREADME, encoding="utf-8").read()
if SENT in cur:
    print("[closeout] %s already carries the CERTIFIED row — skipped" % rel(CREADME))
else:
    outlines = []
    hits = 0
    for line in cur.splitlines():
        if line.startswith(row_prefix):
            hits += 1
            if ("LIVE" not in line) and ("ALL 8 SLICES COMPLETE" not in line):
                die("campaign-I row is in an unrecognised state: %r" % line[:80])
            outlines.append(new_row)
        else:
            outlines.append(line)
    if hits != 1:
        die("expected exactly one campaign-I row in campaigns/README.md, found %d" % hits)
    write(CREADME, "\n".join(outlines) + ("\n" if cur.endswith("\n") else ""))
    print("[closeout] %s — LIVE row flipped to CERTIFIED" % rel(CREADME))

# ========================================================================
# 4c. uc/README.md — headline: the certificate exists at t = psi + 1e-4
# =======================================================================
block = []
block.append("> **CERTIFIED (%s) — a machine-checked certificate exists at the" % DATE)
block.append("> certified rational t = %s ≥ ψ + 10⁻⁴.** Campaign" % launch["t_decimal"])
block.append("> [`cert3_%s_%s`](campaigns/cert3_%s_%s/) is closed: its frozen"
             % (launch["campaign_id"], launch["code_sha256"][:12],
                launch["campaign_id"], launch["code_sha256"][:12]))
block.append("> collector replayed all eight committed proof traces and re-proved")
block.append("> every discharge step, then printed `COMPOSITE CERTIFICATE`:")
block.append("> **Φ ≥ 0 on the feasible 5-parameter family at certified rational")
block.append("> t = %s ≥ exact ψ + %s over w ∈ [1/2,1]**"
             % (launch["t_decimal"], launch["level_offset"]))
block.append("> (the proved orbit-swap symmetry covers w ∈ [0,1]); %s boxes with" % fmt(total))
block.append("> zero residual leaves. Collector (`cert3_collect.py`) SHA-256")
block.append("> `%s`. Verbatim output:" % checker)
block.append("> [`campaigns/CERTIFICATE.txt`](campaigns/CERTIFICATE.txt). %s" % SENT)
block.append("")

cur = open(UREADME, encoding="utf-8").read()
if SENT in cur:
    print("[closeout] %s headline already updated — skipped" % rel(UREADME))
else:
    lines = cur.splitlines(keepends=True)
    h1 = None
    for i, line in enumerate(lines):
        if line.startswith("# "):
            h1 = i
            break
    if h1 is None:
        die("uc/README.md has no level-1 heading")
    j = h1 + 1
    if j < len(lines) and lines[j].strip() == "":
        j += 1
    lines[j:j] = [b + "\n" for b in block]
    write(UREADME, "".join(lines))
    print("[closeout] %s — certificate headline inserted" % rel(UREADME))

print("[closeout] ledger + README updates complete")
PY_MUTATE

say "closeout complete: certificate recorded in PROGRESS.md, campaigns/README.md, uc/README.md; CERTIFICATE.txt archived"
say "manual remainder (see uc/CERT_CLOSEOUT.md): PROOF.md/paper [I-FINAL] placeholders, ANNOUNCEMENT.md, ladder row"
exit 0
