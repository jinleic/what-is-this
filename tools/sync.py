#!/usr/bin/env python3
"""what-is-this: passive progress tracker for the math/ campaigns.

Crawls ~/jinleic-workspace/math, extracts per-project "latest progress"
(project README head + newest matching block from math/PROGRESS.md),
copies a curated set of small core artifacts (papers, conclusions, core
code) into materials/<slug>/, and regenerates the GitHub Pages site
(index.html, problems/*.html, README.md, session-map.md).

Size discipline: only the curated files below are ever copied; anything
missing on disk is skipped; a hard per-file cap protects the repo.

Usage (from a checkout of this repo on the workspace machine):
    python3 tools/sync.py [--limit-kb N]
Then: git add -A && git commit && git push
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------- topology
TOOLS = Path(__file__).resolve()
REPO = TOOLS.parents[1]                      # repos/what-is-this
WORKSPACE = TOOLS.parents[3]                 # ~/jinleic-workspace
MATH = WORKSPACE / "math"
PROGRESS = MATH / "PROGRESS.md"
LIMIT_KB = 5120                              # per-file cap (bytes)

# ---------------------------------------------------------------- projects
# slug -> metadata. sessions: omp session ids that own this project.
PROJECTS = [
    dict(slug="uc", name="Union-closed sets conjecture (Frankl)",
         dirname="uc", session="019ff0a0-6557-7000-8b8c-9bdbccb719e4",
         keywords=("cert3", "union-closed", "Cambie", "Liu", "Frankl", "pscil"),
         anchor="Liu.s Hypothesis 1 is proved",
         papers=("README.md", "ANNOUNCEMENT.md", "PROOF.md",
                 "paper/main.pdf", "paper/main.tex", "paper/refs.bib")),
    dict(slug="ising3d", name="Three-dimensional Ising model (exact solution)",
         dirname="ising3d", session="019ff186-f657-7000-85d0-6d52723bfb5c",
         keywords=("ising3d", "Ising", "wave 1", "Gaussian", "Peierls"),
         anchor="ising3d waves 10",
         papers=("README.md", "research_log.md", "problem_specification.md",
                 "reports/final_technical_report.md",
                 "deliverables/wave18_2026-08-21/papers/main_exact_3d_ising.pdf",
                 "deliverables/wave17_2026-08-20/papers/paperA_no_go_mechanisms.pdf",
                 "deliverables/wave17_2026-08-20/papers/paperB_series_bounds.pdf",
                 "deliverables/wave17_2026-08-20/papers/paperC_layer_DLA.pdf")),
    dict(slug="qec", name="Quantum LDPC codes: exact distance & co-design",
         dirname="qec", session="019ff1a6-af1a-7000-b0e9-223c45ae0a3c",
         keywords=("qec", "QEC", "LDPC", "syzygy", "circuit distance",
                   "logical qubits", "EXP-0", "dressing"),
         anchor=None,
         papers=("README.md",
                 "reports/paper_pbb_nogo.pdf",
                 "reports/fig_trade_law.pdf",
                 "reports/fig_rate_distance.pdf")),
    dict(slug="kobon", name="Kobon triangle problem",
         dirname="kobon", session="019ff95f-d695-7000-b3d0-36dd86537081",
         keywords=("kobon", "KOBON", "K_gen", "triangle", "drat"),
         anchor="KOBON (2026-08-21/22)",
         papers=("README.md", "report.md", "paper_capacity.md",
                 "paper_kobon_2026-08.md", "AUDIT.md",
                 "paper/kobon_broad_capacity.pdf", "paper/kobon_broad_capacity.tex",
                 "engine.py", "release/kobon-2026-08.zip")),
    dict(slug="h10q", name="Hilbert's tenth problem over Q",
         dirname="h10q", session="019ff957-1ed4-7000-8941-f92f3002c676",
         keywords=("h10q", "Hilbert", "L18", "L13", "Schinzel"),
         anchor="L18 step-(ii)",
         papers=("README.md", "RESULTS.md", "THEOREMS.md", "CONDITIONAL.md",
                 "deliverables/2026-08-19/papers/main-conditional-forall6.pdf",
                 "deliverables/2026-08-19/papers/main-conditional-forall6.tex",
                 "deliverables/2026-08-19/papers/companion-verification.pdf",
                 "deliverables/2026-08-19/papers/companion-verification.tex")),
    dict(slug="r55", name="Ramsey number R(5,5)",
         dirname="r55", session="019ffb23-9512-7000-8c79-32667475fc6c",
         keywords=("r55", "R(5,5)", "R(4,5)", "Ramsey", "census", "stratum"),
         anchor=None,
         papers=("README.md", "paper/README.md", "paper/VERIFY.md",
                 "paper/engstrom_mixed_n45.pdf", "paper/engstrom_mixed_n45.tex")),
    dict(slug="zeta5", name="Irrationality of zeta(5) (Apéry-style)",
         dirname="zeta5", session="019ffb23-9512-7000-8c79-32667475fc6c",
         keywords=("zeta5", "Zudilin", "Apery", "Apéry", "zeta(5)", "recsearch"),
         anchor=None,
         papers=("README.md",)),
    dict(slug="ns", name="Navier-Stokes existence & smoothness (route work)",
         dirname="ns", session="019ff0a0-6557-7000-8b8c-9bdbccb719e4 (prelude)",
         keywords=("Navier", "Stokes", "DSS", "self-similar", "Chae"),
         anchor=None,
         papers=("README.md",)),
    dict(slug="ccf", name="Cordoba-Cordoba-Fontelos 1D inviscid model",
         dirname="ccf", session="019ff0a0-6557-7000-8b8c-9bdbccb719e4 (subagents)",
         keywords=("ccf", "Cordoba", "Córdoba", "tail_basis"),
         anchor="tail_basis",
         papers=("README.md",)),
]

HEADER_RE = re.compile(r"^(#{1,3})\s")


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def read_head(path: Path, max_lines: int = 100, max_chars: int = 8000) -> str:
    if not path.exists():
        return ""
    text = path.read_text(errors="replace")
    lines = text.splitlines()[:max_lines]
    out = "\n".join(lines)
    return out[:max_chars]


def readme_tagline(path: Path) -> str:
    """First non-empty, non-heading line of a README == one-line description."""
    if not path.exists():
        return ""
    for line in path.read_text(errors="replace").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            return s[:240]
    return ""


# ---------------------------------------------------------------- progress
def newest_progress_block(project: dict) -> tuple[str, str]:
    """Newest PROGRESS.md block for a project.

    Precedence (PROGRESS.md is append-order newest-first, but campaigns
    interleave, so pure keyword-first-match is unreliable):
      1. block containing the project's anchor regex (if set)
      2. block whose header line contains the slug
      3. block whose body contains the slug as a path ("<slug>/")
      4. first block matching the keyword list
    """
    if not PROGRESS.exists():
        return "", ""
    lines = PROGRESS.read_text(errors="replace").splitlines()
    blocks: list[tuple[int, str, list[str]]] = []
    start = header = None
    body: list[str] = []
    for i, ln in enumerate(lines):
        m = HEADER_RE.match(ln)
        if m and not (start is None and m.group(1) == "#"):
            if start is not None:
                blocks.append((start, header, body))
            start, header, body = i, ln, []
        elif start is not None:
            body.append(ln)
    if start is not None:
        blocks.append((start, header, body))

    slug = project["slug"]
    anchor = project.get("anchor")
    for _, hdr, body in blocks:
        if not anchor:
            continue
        text = "\n".join([hdr] + body)
        if len(text) < 300:
            continue
        if re.search(anchor, text, re.I):
            return hdr, text[:6500]
    for _, hdr, body in blocks:
        h = hdr.lower()
        if (len("\n".join([hdr] + body)) >= 300
                and re.search(rf"(?<![a-z0-9]){slug}(?![a-z0-9])", h)):
            return hdr, "\n".join([hdr] + body)[:6500]
    for _, hdr, body in blocks:
        text = "\n".join([hdr] + body)
        if len(text) < 300:
            continue
        if f"{slug}/" in text.lower():
            return hdr, text[:6500]
    for _, hdr, body in blocks:
        text = "\n".join([hdr] + body)
        if len(text) < 300:
            continue
        if any(k.lower() in text.lower() for k in project["keywords"]):
            return hdr, text[:6500]
    return "", ""


def first_line_nonempty(block: str) -> str:
    for ln in block.splitlines():
        if ln.strip():
            return ln.strip()
    return ""


# ---------------------------------------------------------------- copying
def copy_papers(project: dict) -> list[dict]:
    src = MATH / project["dirname"]
    dest = REPO / "materials" / project["slug"]
    dest.mkdir(parents=True, exist_ok=True)
    copied: list[dict] = []
    for rel in project["papers"]:
        f = src / rel
        if not f.is_file():
            continue
        if f.stat().st_size > LIMIT_KB * 1024:
            print(f"  skip (>{LIMIT_KB}KB): {project['slug']}/{rel} "
                  f"({f.stat().st_size//1024}KB)")
            continue
        out = dest / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(f.read_bytes())
        copied.append(dict(path=rel, size=f.stat().st_size,
                           sha256=hashlib.sha256(f.read_bytes()).hexdigest()[:12]))
    return copied


# ---------------------------------------------------------------- rendering
PAGE_CSS = """
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;
max-width:900px;margin:2rem auto;padding:0 1rem;line-height:1.55;color:#1a1a1a}
h1,h2,h3{line-height:1.25}
pre{background:#f6f8fa;border:1px solid #e1e4e8;border-radius:6px;padding:12px;
overflow-x:auto;font-size:13px;line-height:1.45}
table{border-collapse:collapse;width:100%;margin:1rem 0}
th,td{border:1px solid #d0d7de;padding:6px 10px;text-align:left;font-size:14px}
th{background:#f6f8fa}
a{color:#0969da;text-decoration:none}a:hover{text-decoration:underline}
.tag{font-size:12px;color:#57606a;margin:-0.5rem 0 1.5rem}
.badge{display:inline-block;background:#ddf4ff;border:1px solid #54aeff66;
border-radius:20px;padding:1px 10px;font-size:12px;color:#0969da;margin-right:6px}
.tip{background:#fff8c5;border:1px solid #d4a72c66;border-radius:6px;padding:8px 12px;font-size:13px}
"""


def page(slug: str, title: str, body: str) -> str:
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)} — what-is-this</title>
<style>{PAGE_CSS}</style></head>
<body>
<p><a href="../index.html">&larr; index</a></p>
{body}
<hr><p class="tag">Auto-generated by tools/sync.py · sources under ~/jinleic-workspace/math</p>
</body></html>"""


def build_index(entries: list[dict]) -> str:
    session_map = (REPO / "session-map.md").read_text(errors="replace")
    rows = []
    for e in entries:
        rows.append(
            f"<tr><td><a href=\"{e['page']}\">{esc(e['name'])}</a></td>"
            f"<td><code>{esc(e['dirname'])}</code></td>"
            f"<td><code>{esc(e['session'])}</code></td>"
            f"<td>{esc(e['headline'])}</td>"
            f"<td>{len(e['files'])}</td></tr>")
    table = "\n".join(rows)
    body = f"""
<h1>what-is-this — math campaign progress tracker</h1>
<p class="tag">Public mirror of the research campaigns under
<code>~/jinleic-workspace/math</code>. One page per problem; each page quotes
the project's own README and the newest <code>math/PROGRESS.md</code> entry
verbatim, plus curated papers/artifacts. Regenerated by
<code>tools/sync.py</code> &mdash; a passive tracker, not a mirror of the
workspace.</p>

<table>
<tr><th>Problem</th><th>dir</th><th>session</th><th>Latest headline</th><th>artifacts</th></tr>
{table}
</table>

<h2>How to update</h2>
<pre>cd ~/jinleic-workspace/repos/what-is-this
python3 tools/sync.py        # re-extract progress + copy curated files
git add -A && git commit -m "sync $(date -u +%F)"
git push</pre>

<h2>Session &rarr; project map</h2>
<pre>{esc(session_map)}</pre>
"""
    return body


def build_problem(e: dict) -> str:
    mat = f"../materials/{e['slug']}/"
    fl = "".join(
        f"<li><a href=\"{mat}{esc(f['path'])}\">{esc(f['path'])}</a>"
        f" ({f['size']//1024}KB, sha256 {f['sha256']})</li>"
        for f in e["files"])
    if e["progress_header"]:
        prog = f"<h2>Latest progress — verbatim from math/PROGRESS.md</h2>" \
               f"<p class=\"tag\">Newest entry matching this project</p>" \
               f"<pre>{esc(e['progress'])}</pre>"
    else:
        prog = "<p class=\"tip\">No PROGRESS.md entry matched this project's keywords (minor project — README shown instead).</p>"
    readme = f"<h2>Project README (verbatim, first lines)</h2><pre>{esc(e['readme_head'])}</pre>" \
        if e["readme_head"] else ""
    files = f"<h2>Curated artifacts</h2><ul>{fl}</ul>" if fl else \
        "<p class=\"tip\">No curated artifacts copied yet (project may have none).</p>"
    return f"""
<h1>{esc(e['name'])}</h1>
<p class="tag">source: <code>math/{esc(e['dirname'])}</code> ·
session: <code>{esc(e['session'])}</code> ·
synced {esc(e['synced'])}</p>
{prog}
{readme}
{files}
"""


def build_readme(entries: list[dict]) -> str:
    rows = [f"| {e['name']} | `{e['dirname']}` | `{e['session']}` | {e['headline']} |"
            for e in entries]
    return f"""# what-is-this — math campaign progress tracker

Public progress tracker for the research campaigns under `~/jinleic-workspace/math`.
Regenerated by `tools/sync.py` (passive: quotes each project's own README and the
newest `math/PROGRESS.md` entry verbatim; copies only curated core artifacts).

Live site: https://jinleic.github.io/what-is-this/

## Problems

| Problem | dir | session | Latest headline |
|---|---|---|---|
{chr(10).join(rows)}

## Update

```sh
python3 tools/sync.py
git add -A && git commit -m "sync $(date -u +%F)"
git push
```

See `session-map.md` for the session → project map and `tools/sync.py` for the
curated file list (size-capped; big campaign data lives in the workspace only).
"""


def build_session_map() -> str:
    return """# Session → project map

omp session IDs in `~/.omp/agent/sessions/-jinleic-workspace/<timestamp>_<id>.jsonl`
(head record: title + cwd; first user message names the problem).

| Session | Started | Problem |
|---|---|---|
| `019ff0a0-6557-7000-8b8c-9bdbccb719e4` | 2026-08-11 | union-closed (uc); Millennium/NS prelude; ccf subagents |
| `019ff186-f657-7000-85d0-6d52723bfb5c` | 2026-08-11 | 3D Ising (ising3d) |
| `019ff1a6-af1a-7000-b0e9-223c45ae0a3c` | 2026-08-11 | quantum LDPC / QEC |
| `019ff95f-d695-7000-b3d0-36dd86537081` | 2026-08-13 | Kobon triangles |
| `019ff957-1ed4-7000-8941-f92f3002c676` | 2026-08-13 | Hilbert's tenth over Q (h10q) |
| `019ffb23-9512-7000-8c79-32667475fc6c` | 2026-08-13 | R(5,5) + zeta(5) triage (r55, zeta5) |

Note: `ns` and `ccf` have no dedicated root session (ns was the opening scan of
`019ff0a0`; ccf ran as subagents inside it). Restored 2026-08-15 from the
migration bundle; audit: `SessionIndex` — `019ff0a0` is a union-closed cluster.
"""


# ---------------------------------------------------------------- main
def main() -> int:
    global LIMIT_KB
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit-kb", type=int, default=LIMIT_KB // 1024)
    args = ap.parse_args()
    LIMIT_KB = args.limit_kb * 1024

    if not MATH.is_dir():
        print(f"source math/ not found at {MATH}", file=sys.stderr)
        return 1

    (REPO / "problems").mkdir(exist_ok=True)
    (REPO / "materials").mkdir(exist_ok=True)
    (REPO / ".nojekyll").write_text("")

    ts = now()
    entries = []
    total_bytes = 0
    for p in PROJECTS:
        slug = p["slug"]
        src = MATH / p["dirname"]
        readme_path = src / "README.md"
        readme_head = read_head(readme_path)
        tagline = readme_tagline(readme_path)
        header, block = newest_progress_block(p)
        files = copy_papers(p)
        total_bytes += sum(f["size"] for f in files)
        headline = first_line_nonempty(block)
        # prefer the anchor block name when present; fall back to first line
        for probe in (p.get("anchor"),):
            if probe and block and re.search(probe, block, re.I):
                m = re.search(r"(?m)^#{1,3}\s+.*$", block)
                headline = (m.group(0) if m else first_line_nonempty(block))[:180]
                break
        e = dict(slug=slug, name=p["name"], dirname=p["dirname"],
                 session=p["session"], tagline=tagline, headline=headline[:180],
                 progress_header=header, progress=block, readme_head=readme_head,
                 files=files, synced=ts,
                 page=f"problems/{slug}.html")
        entries.append(e)
        (REPO / "problems" / f"{slug}.html").write_text(
            page(slug, p["name"], build_problem(e)))
        print(f"{slug:8s} headline: {headline[:110] or '(none)'}")
        print(f"{'':8s} artifacts: {len(files)} files")

    entries_sorted = sorted(entries, key=lambda e: e["name"].lower())
    (REPO / "session-map.md").write_text(build_session_map())
    (REPO / "README.md").write_text(build_readme(entries_sorted))
    (REPO / "index.html").write_text(
        page("index", "what-is-this — math campaign progress tracker",
             build_index(entries_sorted)))

    manifest = dict(generated=ts,
                    limit_kb=args.limit_kb,
                    total_copied_bytes=total_bytes,
                    projects={e["slug"]: dict(
                        synced=e["synced"], progress_header=e["progress_header"],
                        readme_head_lines=len(e["readme_head"].splitlines()),
                        copied_bytes=sum(f["size"] for f in e["files"]),
                        files=e["files"],
                        headline=e["headline"]) for e in entries})
    (REPO / "sync.json").write_text(json.dumps(manifest, indent=1))

    print(f"\nsite regenerated at {ts}")
    print(f"total curated bytes: {total_bytes//1024}KB (~{total_bytes/1024/1024:.1f}MB)")
    print("commit + push to publish; see README.md for the update loop")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())