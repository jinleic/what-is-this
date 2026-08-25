#!/usr/bin/env python3
"""what-is-this: passive progress tracker + static progress mirror.

Mirrors the math campaign workspace into the repo as a public progress
mirror (no website hosting):

- one page per problem, quoting the project's README head and the newest
  matching ``math/PROGRESS.md`` entry; both are rendered client-side
  (markdown via marked.js, LaTeX via MathJax with math-span protection);
- mirrors each project's core code and small artifacts (per-file size cap,
  very large result dirs excluded);
- incremental: a file is rewritten only when its content changed, and each
  project keeps a "last change" timestamp that survives no-op runs, so a
  publish is always a minimal diff;
- private-keyword gate: words listed in a PRIVATE file that lives OUTSIDE
  this repository must never appear in anything published.  Source files
  containing them are withheld, rendered text is redacted, and the run
  fails (non-zero exit) if any hit survives in the output tree.  The run
  also refuses to start when the private list is missing (fail-closed).

Usage: python3 tools/sync.py [--limit-kb N]
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------- topology
TOOLS = Path(__file__).resolve()
REPO = TOOLS.parents[1]                      # repos/what-is-this
WORKSPACE = TOOLS.parents[3]                 # ~/jinleic-workspace
MATH = WORKSPACE / "math"
PROGRESS = MATH / "PROGRESS.md"
LIMIT_KB = 5120                              # per-file cap (KB)
PRIVATE_WORDS = WORKSPACE / ".sync-banned.txt"   # NEVER inside the repo

# ---------------------------------------------------------------- projects
PROJECTS = [
    dict(slug="uc", name="Union-closed sets conjecture (Frankl)",
         dirname="uc", session="019ff0a0-6557-7000-8b8c-9bdbccb719e4",
         keywords=("cert3", "union-closed", "Cambie", "Liu", "Frankl", "pscil"),
         anchor="campaign J:",
         exclude_dirs=("campaigns",),
         papers=("README.md", "ANNOUNCEMENT.md", "PROOF.md",
                 "paper/main.pdf", "paper/main.tex", "paper/refs.bib")),
    dict(slug="ising3d", name="Three-dimensional Ising model (exact solution)",
         dirname="ising3d", session="019ff186-f657-7000-85d0-6d52723bfb5c",
         keywords=("ising3d", "Ising", "wave 1", "Gaussian", "Peierls"),
         anchor="ISING3D WAVE 22",
         exclude_dirs=("results", ".venv"),
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
         exclude_dirs=(".venv", "third_party", "dist"),
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
         anchor="H10/Q",
         exclude_dirs=("data",),
         papers=("README.md", "RESULTS.md", "THEOREMS.md", "CONDITIONAL.md",
                 "deliverables/2026-08-19/papers/main-conditional-forall6.pdf",
                 "deliverables/2026-08-19/papers/main-conditional-forall6.tex",
                 "deliverables/2026-08-19/papers/companion-verification.pdf",
                 "deliverables/2026-08-19/papers/companion-verification.tex")),
    dict(slug="r55", name="Ramsey number R(5,5)",
         dirname="r55", session="019ffb23-9512-7000-8c79-32667475fc6c",
         keywords=("r55", "R(5,5)", "R(4,5)", "Ramsey", "census", "stratum"),
         anchor=None,
         exclude_dirs=("data", "outD"),
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


# ---------------------------------------------------------------- helpers
def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def read_head(path: Path, max_lines: int = 100, max_chars: int = 8000) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(errors="replace").splitlines()[:max_lines]
    return "\n".join(lines)[:max_chars]


def readme_tagline(path: Path) -> str:
    if not path.exists():
        return ""
    for line in path.read_text(errors="replace").splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            return s[:240]
    return ""


def first_line_nonempty(block: str) -> str:
    for ln in block.splitlines():
        if ln.strip():
            return ln.strip()
    return ""


def strip_md(s: str) -> str:
    """Inline-markdown -> plain text (for table cells; keeps $math$)."""
    s = re.sub(r"^#{1,6}\s*", "", s)
    s = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", s)
    return s.replace("**", "").replace("`", "").strip()


def write_if_changed(path: Path, text: str) -> bool:
    data = text.encode()
    if path.exists() and path.read_bytes() == data:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return True


# ---------------------------------------------------------------- privacy
def load_private_words() -> list[str]:
    """Words that must never be published.  Fail-closed."""
    words: list[str] = []
    if PRIVATE_WORDS.exists():
        words = [w.strip().lower()
                 for w in PRIVATE_WORDS.read_text().splitlines() if w.strip()]
    if not words:
        print(f"FATAL: private keyword list missing or empty: {PRIVATE_WORDS}")
        print("refusing to sync (fail-closed); create the file, one word per line.")
        raise SystemExit(2)
    return words


def private_hit(data: bytes, words: list[str]) -> bool:
    low = data.lower()
    return any(w.encode() in low for w in words)


def redact(text: str, words: list[str]) -> str:
    for w in words:
        text = re.sub(re.escape(w), "[private]", text, flags=re.I)
    return text


def final_gate(words: list[str]) -> int:
    """Scan everything that would be published; non-zero exit on any hit."""
    bad = []
    for f in REPO.rglob("*"):
        if not f.is_file() or ".git" in f.parts:
            continue
        if private_hit(f.read_bytes(), words):
            bad.append(str(f.relative_to(REPO)))
    if bad:
        print("PRIVATE-KEYWORD GATE FAILED — do not publish; offending files:")
        for b in bad:
            print("   ", b)
        return 1
    print("private-keyword gate: clean")
    return 0


# ---------------------------------------------------------------- progress
def newest_progress_block(project: dict) -> tuple[str, str]:
    """Newest PROGRESS.md block for a project.

    Precedence (PROGRESS.md is append-order newest-first, but campaigns
    interleave, so pure keyword-first-match is unreliable):
      1. block containing the project's anchor regex (if set)
      2. block whose header line contains the slug (word-boundary)
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
    if anchor:
        for _, hdr, bd in blocks:
            text = "\n".join([hdr] + bd)
            if len(text) >= 300 and re.search(anchor, text, re.I):
                return hdr, text[:6500]
    for _, hdr, bd in blocks:
        text = "\n".join([hdr] + bd)
        if (len(text) >= 300
                and re.search(rf"(?<![a-z0-9]){slug}(?![a-z0-9])", hdr.lower())):
            return hdr, text[:6500]
    for _, hdr, bd in blocks:
        text = "\n".join([hdr] + bd)
        if len(text) >= 300 and f"{slug}/" in text.lower():
            return hdr, text[:6500]
    for _, hdr, bd in blocks:
        text = "\n".join([hdr] + bd)
        if len(text) >= 300 and any(k.lower() in text.lower()
                                    for k in project["keywords"]):
            return hdr, text[:6500]
    return "", ""


# ---------------------------------------------------------------- mirror
JUNK_DIRS = {"__pycache__", ".venv", "venv", ".git", "node_modules",
             ".ipynb_checkpoints", ".pytest_cache", ".mypy_cache",
             ".ruff_cache", "build", "dist", "egg-info"}
JUNK_EXTS = {".pyc", ".pyo"}


def collect_sources(project: dict) -> dict[str, Path]:
    """rel-path -> source file.  Curated papers first, then the tree walk."""
    src = MATH / project["dirname"]
    out: dict[str, Path] = {}
    for rel in project["papers"]:
        f = src / rel
        if not f.is_file():
            continue
        if f.stat().st_size > LIMIT_KB * 1024:
            print(f"  skip (>{LIMIT_KB}KB): {project['slug']}/{rel}")
            continue
        out[rel] = f
    exclude = set(project.get("exclude_dirs", ()))
    if not src.is_dir():
        return out
    for f in sorted(src.rglob("*")):
        if not f.is_file():
            continue
        rel = f.relative_to(src)
        if any(p in JUNK_DIRS or p in exclude for p in rel.parts[:-1]):
            continue
        if f.suffix in JUNK_EXTS or f.name == ".DS_Store":
            continue
        if f.stat().st_size > LIMIT_KB * 1024:
            continue
        out.setdefault(str(rel), f)
    return out


def sync_materials(project: dict, words: list[str]) -> tuple[list[dict], dict]:
    """Incrementally mirror one project; returns (files, stats)."""
    dest = REPO / "materials" / project["slug"]
    dest.mkdir(parents=True, exist_ok=True)
    files: list[dict] = []
    stats = dict(written=0, unchanged=0, pruned=0, private=0)
    kept: set[str] = set()
    for rel, srcf in collect_sources(project).items():
        data = srcf.read_bytes()
        if private_hit(data, words):
            stats["private"] += 1
            continue
        sha = hashlib.sha256(data).hexdigest()[:12]
        out = dest / rel
        kept.add(rel)
        if (out.exists() and out.stat().st_size == len(data)
                and hashlib.sha256(out.read_bytes()).hexdigest()[:12] == sha):
            stats["unchanged"] += 1
        else:
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(data)
            stats["written"] += 1
        files.append(dict(path=rel, size=len(data), sha256=sha))
    # prune anything that no longer belongs to the mirror set
    for f in sorted(dest.rglob("*")):
        if f.is_file() and str(f.relative_to(dest)) not in kept:
            f.unlink()
            stats["pruned"] += 1
    for d in sorted((p for p in dest.rglob("*") if p.is_dir()), reverse=True):
        try:
            d.rmdir()
        except OSError:
            pass
    files.sort(key=lambda x: x["path"])
    return files, stats


# ---------------------------------------------------------------- rendering
PAGE_CSS = """
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;
max-width:960px;margin:2rem auto;padding:0 1rem;line-height:1.55;color:#1f2328}
h1,h2,h3{line-height:1.25}
pre{background:#f6f8fa;border:1px solid #e1e4e8;border-radius:6px;padding:12px;
overflow-x:auto;font-size:13px;line-height:1.45}
table{border-collapse:collapse;width:100%;margin:1rem 0}
th,td{border:1px solid #d0d7de;padding:6px 10px;text-align:left;font-size:14px;vertical-align:top}
th{background:#f6f8fa}
a{color:#0969da;text-decoration:none}a:hover{text-decoration:underline}
.tag{font-size:12px;color:#57606a;margin:-0.4rem 0 1rem}
.tip{background:#fff8c5;border:1px solid #d4a72c66;border-radius:6px;padding:8px 12px;font-size:13px}
.markdown-body{box-sizing:border-box;background:transparent;font-size:15px;margin:0 0 1.2rem}
details{border:1px solid #d0d7de;border-radius:6px;padding:6px 12px;margin:8px 0;background:#fafbfc}
summary{cursor:pointer;font-size:14px}
ul.files{margin:8px 0;padding-left:20px;font-size:13px}
ul.files li{margin:2px 0}
.fsize{color:#57606a;font-size:11px}
mjx-container{overflow-x:auto;max-width:100%}
"""

HEAD_STATIC = """
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/github-markdown-css@5.5.1/github-markdown-light.min.css">
<script>
window.MathJax = {
  tex: {
    inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
    displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
    processEscapes: true,
    packages: {'[+]': ['color', 'ams']}
  }
};
</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js"></script>
<script src="https://cdn.jsdelivr.net/npm/marked@12.0.2/marked.min.js"></script>
<script>
(function () {
  function renderMD(el, md) {
    var math = [];
    function stash(re) {
      md = md.replace(re, function (m) {
        math.push(m);
        return '@@MATH' + (math.length - 1) + '@@';
      });
    }
    stash(/\\$\\$[\\s\\S]+?\\$\\$/g);
    stash(/\\\\\\[[\\s\\S]+?\\\\\\]/g);
    stash(/\\\\\\((?:[^\\\\]|\\\\[^)])+?\\\\\\)/g);
    stash(/\\$(?=\\S)(?:[^$\\n]*?\\S)?\\$/g);
    var out = marked.parse(md, { mangle: false, headerIds: false });
    out = out.replace(/@@MATH(\\d+)@@/g, function (_, i) {
      return math[+i].replace(/&/g, '&amp;')
                     .replace(/</g, '&lt;').replace(/>/g, '&gt;');
    });
    el.innerHTML = out;
  }
  document.addEventListener('DOMContentLoaded', function () {
    function rewriteLinks(el) {
      var base = (window.__LINKS__ || {})[el.id];
      if (!base) return;
      var nodes = el.querySelectorAll('a[href], img[src]');
      for (var i = 0; i < nodes.length; i++) {
        var n = nodes[i];
        var attr = n.tagName === 'IMG' ? 'src' : 'href';
        var h = n.getAttribute(attr);
        if (!h) continue;
        if (/^(?:https?:|mailto:|tel:|#|\/\/|\/)/i.test(h)) continue;
        n.setAttribute(attr, base + h);
      }
    }
    var blobs = window.__MD__ || {};
    Object.keys(blobs).forEach(function (id) {
      var el = document.getElementById(id);
      if (!el) return;
      try { renderMD(el, blobs[id]); rewriteLinks(el); }
      catch (err) {
        var p = document.createElement('pre');
        p.textContent = blobs[id];
        el.innerHTML = '';
        el.appendChild(p);
      }
    });
    if (window.MathJax) {
      if (MathJax.typesetPromise) { MathJax.typesetPromise(); }
      else if (MathJax.startup && MathJax.startup.promise) {
        MathJax.startup.promise.then(function () { MathJax.typesetPromise(); });
      }
    }
  });
})();
</script>
"""


def page(title: str, body: str, md_map: dict[str, str], depth: int = 0,
         links: dict[str, str] | None = None) -> str:
    md_json = json.dumps(md_map, ensure_ascii=False).replace("<", "\\u003c")
    lk_json = json.dumps(links or {}, ensure_ascii=False)
    home = "../index.html" if depth else "index.html"
    back = f'<p><a href="{home}">&larr; index</a></p>' if depth else ""
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)} — what-is-this</title>
<style>{PAGE_CSS}</style>{HEAD_STATIC}</head>
<body>
{back}
{body}
<script>window.__MD__ = {md_json};</script>
<script>window.__LINKS__ = {lk_json};</script>
<hr><p class="tag">Auto-generated by tools/sync.py · sources under ~/jinleic-workspace/math</p>
</body></html>"""


def build_problem(e: dict) -> tuple[str, dict[str, str]]:
    mat = f"../materials/{e['slug']}/"
    groups: dict[str, list[dict]] = {}
    for f in e["files"]:
        top = f["path"].split("/", 1)[0] if "/" in f["path"] else "(root)"
        groups.setdefault(top, []).append(f)
    parts = []
    for g in sorted(groups, key=lambda k: (k != "(root)", k.lower())):
        fs = groups[g]
        mb = sum(x["size"] for x in fs) / 1e6
        lis = "".join(
            f'<li><a href="{mat}{esc(x["path"])}">{esc(x["path"])}</a>'
            f' <span class="fsize">{max(x["size"] // 1024, 1)} KB</span></li>'
            for x in fs)
        label = "project root" if g == "(root)" else g + "/"
        op = " open" if g == "(root)" else ""
        parts.append(f'<details{op}><summary><code>{esc(label)}</code>'
                     f' — {len(fs)} files, {mb:.1f} MB</summary>'
                     f'<ul class="files">{lis}</ul></details>')
    files_html = "\n".join(parts) or '<p class="tip">No artifacts copied.</p>'

    md_map: dict[str, str] = {}
    if e["progress"]:
        md_map["md-progress"] = e["progress"]
        prog = ('<h2>Latest progress</h2><p class="tag">newest matching entry, '
                'quoted from <code>math/PROGRESS.md</code></p>'
                '<div class="markdown-body" id="md-progress"></div>')
    else:
        prog = '<p class="tip">No PROGRESS.md entry matched this project.</p>'
    if e["readme_head"]:
        md_map["md-readme"] = e["readme_head"]
        readme = ('<h2>Project README</h2><p class="tag">first lines, verbatim</p>'
                  '<div class="markdown-body" id="md-readme"></div>')
    else:
        readme = ""
    body = f"""<h1>{esc(e['name'])}</h1>
<p class="tag">source <code>math/{esc(e['dirname'])}</code> ·
session <code>{esc(e['session'])}</code> ·
last change {esc(e['updated'])}</p>
{prog}
{readme}
<h2>Code &amp; artifacts</h2>
<p class="tag">mirrored from the project tree; very large result files and campaign data excluded</p>
{files_html}"""
    return body, md_map, {
        "md-progress": "../materials/",      # PROGRESS.md lives at math/ root
        "md-readme": f"../materials/{e['slug']}/",
    }


def build_index(entries: list[dict], session_md: str) -> tuple[str, dict[str, str]]:
    rows = []
    for e in entries:
        rows.append(
            "<tr>"
            f'<td><a href="{e["page"]}">{esc(e["name"])}</a></td>'
            f'<td><code>{esc(e["dirname"])}</code></td>'
            f'<td>{esc(e["updated"][:10])}</td>'
            f'<td>{esc(strip_md(e["headline"]))}</td>'
            f'<td>{len(e["files"])}</td>'
            "</tr>")
    table = "\n".join(rows)
    body = f"""<h1>what-is-this — math campaign progress tracker</h1>
<p class="tag">Public mirror of the research campaigns under
<code>~/jinleic-workspace/math</code>. One page per problem: the project's
README head and the newest matching <code>math/PROGRESS.md</code> entry,
quoted verbatim and rendered (markdown + LaTeX), plus a mirror of the
project's core code and small artifacts. Very large result data stays in the
workspace. Incremental: pages and files are rewritten only when content
changes.</p>

<table>
<tr><th>Problem</th><th>dir</th><th>last change</th><th>latest headline</th><th>files</th></tr>
{table}
</table>

<h2>How to update</h2>
<pre>cd ~/jinleic-workspace/repos/what-is-this
python3 tools/sync.py        # incremental; fails closed on the private-keyword gate
git add -A && git commit -m "sync $(date -u +%F)" && git push</pre>

<h2>Sessions</h2>
<div class="markdown-body" id="md-sessions"></div>"""
    return body, {"md-sessions": session_md}, {"md-sessions": ""}


def build_readme(entries: list[dict]) -> str:
    rows = []
    for e in entries:
        h = strip_md(e["headline"]).replace("|", "\\|")
        rows.append(f"| [{e['name']}](problems/{e['slug']}.html) "
                    f"| `{e['dirname']}` | {e['updated'][:10]} | {h} |")
    nl = "\n"
    return f"""# what-is-this — math campaign progress tracker

Public progress tracker for the research campaigns under `~/jinleic-workspace/math`.
Regenerated by `tools/sync.py` (passive: quotes each project's own README and the
newest `math/PROGRESS.md` entry verbatim; mirrors core code + small artifacts,
never multi-GB campaign data or single files over the size cap; incremental —
only changed files are rewritten; private-keyword gate blocks publishing when
any word from a private, out-of-repo list appears in the output).

## Problems

| Problem | dir | last change | latest headline |
|---|---|---|---|
{nl.join(rows)}

## Update

```sh
python3 tools/sync.py
git add -A && git commit -m "sync $(date -u +%F)"
git push
```

See `session-map.md` for the session → project map and `tools/sync.py` for the
mirror rules (per-file size cap, excluded result dirs, private-keyword gate).
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

`ns` and `ccf` have no dedicated root session: NS route work happened inside the
`019ff0a0` prelude, and ccf inside its subagents.
"""


# ---------------------------------------------------------------- main
def main() -> int:
    global LIMIT_KB
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit-kb", type=int, default=LIMIT_KB,
                    help="per-file size cap in KB")
    args = ap.parse_args()
    LIMIT_KB = args.limit_kb

    words = load_private_words()
    ts = now()
    sj = REPO / "sync.json"
    prev: dict = {}
    if sj.exists():
        try:
            prev = json.loads(sj.read_text())
        except Exception:
            prev = {}
    prev_projects = prev.get("projects", {})

    entries: list[dict] = []
    total_bytes = 0
    pages_rewritten = 0
    for p in PROJECTS:
        slug = p["slug"]
        readme_path = MATH / p["dirname"] / "README.md"
        readme_head = redact(read_head(readme_path), words)
        tagline = redact(readme_tagline(readme_path), words)
        header, block = newest_progress_block(p)
        header = redact(header, words)
        block = redact(block, words)
        files, st = sync_materials(p, words)
        total_bytes += sum(f["size"] for f in files)
        m = re.search(r"(?m)^#{1,3}\s+.*$", block)
        headline = (m.group(0) if m else first_line_nonempty(block))[:180]
        sig = hashlib.sha256(
            (block + "\x00" + readme_head + "\x00"
             + json.dumps(files, sort_keys=True)).encode()).hexdigest()[:16]
        pv = prev_projects.get(slug, {})
        updated = pv.get("updated", ts) if pv.get("sig") == sig else ts
        e = dict(slug=slug, name=p["name"], dirname=p["dirname"],
                 session=p["session"], tagline=tagline, headline=headline,
                 progress_header=header, progress=block,
                 readme_head=readme_head, files=files,
                 updated=updated, sig=sig, page=f"problems/{slug}.html")
        entries.append(e)
        body, md_map, links = build_problem(e)
        if write_if_changed(REPO / "problems" / f"{slug}.html",
                            page(p["name"], body, md_map, depth=1, links=links)):
            pages_rewritten += 1
        print(f"{slug:8s} {headline[:100] or '(no entry)'}")
        print(f"{'':8s} files={len(files)} "
              f"({sum(f['size'] for f in files)/1e6:.1f} MB): "
              f"{st['written']} written, {st['unchanged']} unchanged, "
              f"{st['pruned']} pruned, {st['private']} withheld (private)")

    entries_sorted = sorted(entries, key=lambda e: e["name"].lower())
    session_md = build_session_map()
    write_if_changed(REPO / "session-map.md", session_md)
    write_if_changed(REPO / "README.md", build_readme(entries_sorted))
    idx_body, idx_md, idx_links = build_index(entries_sorted, session_md)
    if write_if_changed(REPO / "index.html",
                        page("what-is-this — math campaign progress tracker",
                             idx_body, idx_md, depth=0, links=idx_links)):
        pages_rewritten += 1

    manifest = dict(
        updated=max((e["updated"] for e in entries), default=ts),
        limit_kb=LIMIT_KB,
        projects={e["slug"]: dict(updated=e["updated"], sig=e["sig"],
                                  headline=e["headline"],
                                  progress_header=e["progress_header"],
                                  bytes=sum(f["size"] for f in e["files"]),
                                  files=e["files"]) for e in entries})
    write_if_changed(sj, json.dumps(manifest, indent=1, sort_keys=True,
                                    ensure_ascii=False))

    print(f"\ntotal mirrored: {total_bytes/1e6:.1f} MB; "
          f"pages rewritten this run: {pages_rewritten}")
    rc = final_gate(words)
    if rc == 0:
        print("clean — sync with: git add -A && git commit && git push")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
