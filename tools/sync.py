#!/usr/bin/env python3
"""what-is-this: passive research progress tracker + static artifact mirror.

Mirrors the workspace's mathematics, physics, theoretical-computer-science,
and quantitative-trading research into this public repository:

- one page per target, quoting the target README head and the newest matching
  domain ``PROGRESS.md`` entry; both render client-side (Markdown + LaTeX);
- one curated quantitative-trading page sourced from allowlisted status fields;
- one copy of each domain's root README/results/progress documents;
- each target's core code and small artifacts, subject to a per-file size cap
  and explicit exclusions for raw, provisional, or third-party material;
- incremental writes and stable per-target "last change" timestamps;
- a fail-closed publication gate: private keywords live outside this repo,
  high-confidence secret shapes are blocked, and unsafe source files are
  withheld before anything can be published.

Existing mathematics URLs remain unchanged.

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
REPO = TOOLS.parents[1]                          # repos/what-is-this
WORKSPACE = TOOLS.parents[3]                     # ~/jinleic-workspace
LIMIT_KB = 5120                                  # per-file cap (KB)
PRIVATE_WORDS = WORKSPACE / ".sync-banned.txt"   # NEVER inside the repo

DOMAINS = {
    "math": dict(
        name="Mathematics",
        root=WORKSPACE / "math",
        progress="PROGRESS.md",
        shared_files=("README.md", "RESULTS.md", "PROGRESS.md",
                      "RESEARCH_STATUS.md"),
    ),
    "physics": dict(
        name="Physics and quantum computing",
        root=WORKSPACE / "physics",
        progress="PROGRESS.md",
        shared_files=("README.md", "RESULTS.md", "PROGRESS.md"),
    ),
    "cs": dict(
        name="Theoretical computer science and information theory",
        root=WORKSPACE / "cs",
        progress="PROGRESS.md",
        shared_files=("README.md", "RESULTS.md", "PROGRESS.md"),
    ),
    "quant-trading": dict(
        name="Quantitative trading research",
        root=WORKSPACE / "quant-trading",
        progress=None,
        shared_files=(),
    ),
}

VERIFIED_TARGET_EXCLUDES = ("scratch",)

# ---------------------------------------------------------------- projects
PROJECTS = [
    dict(slug="e389",
         name="Erdős Problem #389: consecutive-product divisibility",
         dirname="e389", session="01a03c0c-c36d-73b7-96b7-f1203909eb76",
         keywords=("E389", "Erdős Problem #389", "consecutive-product",
                   "compensation-good"),
         anchor=None,
         exclude_dirs=("data",),
         papers=("README.md", "THEOREMS.md",
                 "data/known_witness_verification.json",
                 "data/compensation_run_m27_k50001_h100000000.json",
                 "data/compensation_run_m27_k100050001_h900000000.json",
                 "data/compensation_structure_m1_20_k5000_m27_h200000000.json")),
    dict(slug="liu_h1",
         name="Liu Hypothesis 1: entropy-kernel theorem",
         dirname="LIU_H1", session="01a03b5f-2600-7019-adf6-b03bd015c2ce",
         keywords=("LIU_H1", "Liu Hypothesis 1", "Taylor--Lorentz--Gram"),
         anchor=None,
         exclude_dirs=("literature", "logs"),
         papers=("AUDIT.md", "REPRODUCIBILITY.md",
                 "LITERATURE_ORIGINALITY.md", "paper/main.pdf",
                 "paper/main.tex", "paper/refs.bib",
                 "verification/independent_exact_checker.py",
                 "verification/environment.json")),
    dict(slug="uc", name="Union-closed sets conjecture (Frankl)",
         dirname="uc", session="019ff0a0-6557-7000-8b8c-9bdbccb719e4",
         keywords=("cert3", "union-closed", "Cambie", "Liu", "Frankl", "pscil"),
         anchor="UC-",
         exclude_dirs=("independent-arithmetic",),
         papers=("README.md", "ANNOUNCEMENT.md", "PROOF.md",
                 "paper/main.pdf", "paper/main.tex", "paper/refs.bib")),
    dict(slug="ising3d", name="Three-dimensional Ising model (exact solution)",
         dirname="ising3d", session="019ff186-f657-7000-85d0-6d52723bfb5c",
         keywords=("ising3d", "Ising", "wave 1", "Gaussian", "Peierls"),
         anchor="ISING-W",
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
         exclude_dirs=(".venv", "third_party", "dist", "results", "scratch",
                       "archive"),
         papers=("README.md",
                 "reports/paper_pbb_nogo.pdf",
                 "reports/fig_trade_law.pdf",
                 "reports/fig_rate_distance.pdf")),
    dict(slug="kobon", name="Kobon triangle problem",
         dirname="kobon", session="019ff95f-d695-7000-b3d0-36dd86537081",
         keywords=("kobon", "KOBON", "K_gen", "triangle", "drat"),
         anchor="KOBON",
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
    # Physics targets. Quantum computing is owned here; there is no separate
    # top-level quantum research directory.
    dict(domain="physics", slug="qldpc-dec",
         name="BB quantum-LDPC decoder reproduction",
         dirname="qldpc-dec", session="",
         keywords=("qldpc-dec", "beam search", "GARI", "BP+OSD"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES + ("_gateB_shard_cmd",),
         papers=("README.md",)),
    dict(domain="physics", slug="msd",
         name="Zero-level CCZ reproduction",
         dirname="msd", session="",
         keywords=("msd", "zero-level CCZ", "magic-state"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="physics", slug="shadows",
         name="Contractive-shadow exact census",
         dirname="shadows", session="",
         keywords=("shadows", "contractive-shadow", "Clifford"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="physics", slug="na-compiler",
         name="Certified neutral-atom transport scheduling",
         dirname="na-compiler", session="",
         keywords=("na-compiler", "neutral-atom", "transport scheduling"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="physics", slug="fss-bb",
         name="BB-code finite-size scaling",
         dirname="fss-bb", session="",
         keywords=("fss-bb", "finite-size scaling", "ν anomaly"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="physics", slug="qlops",
         name="Fault-tolerant QLOPS arithmetic audit",
         dirname="qlops", session="",
         keywords=("qlops", "QLOPS", "resource-estimate"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),

    # Theoretical computer science and information theory targets.
    dict(domain="cs", slug="mceliece",
         name='Classic McEliece hold-out "waterfall" dispute',
         dirname="mceliece", session="",
         keywords=("mceliece", "Goppa", "waterfall"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES
                      + ("prim", "preserving_precious"),
         papers=("README.md",)),
    dict(domain="cs", slug="kg",
         name="Grothendieck constant",
         dirname="kg", session="",
         keywords=("kg/", "Grothendieck", "K_G"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="cs", slug="omega",
         name="Matrix multiplication exponent",
         dirname="omega", session="",
         keywords=("omega/", "matrix multiplication exponent", "ω"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="cs", slug="delcap",
         name="Binary deletion-channel capacity",
         dirname="delcap", session="",
         keywords=("delcap", "deletion channel", "Blahut"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="cs", slug="oct-rank",
         name="Real tensor rank of octonion multiplication",
         dirname="oct-rank", session="",
         keywords=("oct-rank", "octonion", "tensor rank"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="cs", slug="rs-pe3d",
         name="Reed–Solomon three-dimensional product expansion",
         dirname="rs-pe3d", session="",
         keywords=("rs-pe3d", "product expansion", "Reed–Solomon"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         papers=("README.md",)),
    dict(domain="cs", slug="mm3",
         name="Additive complexity of rank-23 3×3 matrix multiplication",
         dirname="mm3", session="",
         keywords=("mm3", "rank-23", "55 additions"),
         anchor=None,
         exclude_dirs=VERIFIED_TARGET_EXCLUDES,
         exclude_paths=(".README.md.lock",),
         papers=("README.md",)),

    # Quantitative trading has one curated overview. Raw market data, return
    # series, downloads, provisional outputs, and machine-local paths stay out.
    dict(domain="quant-trading", slug="overview",
         name="Quantitative trading research",
         dirname=".", session="",
         keywords=(), anchor=None, papers=(),
         summary="quant-trading", single_page=True,
         exclude_dirs=("data", "downloads", "output", "outputs", "artifacts",
                       "analysis", "results", "scratch",
                       "decomposition_output"),
         exclude_paths=(
             "workspace.json",
             "explorations/accelerated-intraday-blinding-audit.json",
         ),
         include_exts=(".json", ".md", ".py"),
         withhold_home_paths=True),
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


def project_domain(project: dict) -> str:
    return project.get("domain", "math")


def project_source(project: dict) -> Path:
    return DOMAINS[project_domain(project)]["root"] / project["dirname"]


def project_page(project: dict) -> str:
    if "page" in project:
        return project["page"]
    domain = project_domain(project)
    if project.get("single_page"):
        return f"{domain}/index.html"
    root = "problems" if domain == "math" else domain
    return f"{root}/{project['slug']}.html"


def project_material(project: dict) -> str:
    if "material" in project:
        return project["material"]
    domain = project_domain(project)
    if domain == "math":
        return f"materials/{project['slug']}"
    if project.get("single_page"):
        return f"materials/{domain}"
    return f"materials/{domain}/{project['slug']}"


def project_manifest_key(project: dict) -> str:
    domain = project_domain(project)
    return project["slug"] if domain == "math" else f"{domain}/{project['slug']}"


def project_source_label(project: dict) -> str:
    domain = project_domain(project)
    dirname = project["dirname"]
    return domain if dirname == "." else f"{domain}/{dirname}"


JUNK_DIRS = {"__pycache__", ".venv", "venv", ".git", "node_modules",
             ".ipynb_checkpoints", ".pytest_cache", ".mypy_cache",
             ".ruff_cache", "build", "dist", "egg-info",
             # Campaign run directories are timestamped raw outputs, never
             # curated mirror content (their summaries live in PROGRESS.md).
             "campaigns", "campaigns-smoke"}
# Archives are opaque to the text-based publication gate, so they are never
# walk-mirrored; curated archives ship only via an explicit `papers` entry.
JUNK_EXTS = {".pyc", ".pyo", ".zip", ".tar", ".tgz", ".gz", ".bz2",
             ".xz", ".zst", ".7z", ".rar"}

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


SENSITIVE_PATTERNS = (
    re.compile(rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    re.compile(rb"s[k]-[A-Za-z0-9_-]{20,}"),
    re.compile(rb"g[h][pousr]_[A-Za-z0-9]{20,}"),
    re.compile(rb"A[K]IA[0-9A-Z]{16}"),
    re.compile(rb"https://hooks[.]slack[.]com/services/[A-Za-z0-9/_-]+"),
)
HOME_PATH_RE = re.compile(rb"/(?:Users|home)/[^/\s]+/")


def private_hit(data: bytes, words: list[str]) -> bool:
    low = data.lower()
    return any(w.encode() in low for w in words)


def sensitive_hit(data: bytes) -> bool:
    return any(pattern.search(data) for pattern in SENSITIVE_PATTERNS)


def publication_hit(data: bytes, words: list[str]) -> bool:
    return private_hit(data, words) or sensitive_hit(data)


def source_blocked(project: dict, data: bytes, words: list[str]) -> bool:
    if publication_hit(data, words):
        return True
    return bool(project.get("withhold_home_paths") and HOME_PATH_RE.search(data))


def redact(text: str, words: list[str]) -> str:
    for w in words:
        text = re.sub(re.escape(w), "[private]", text, flags=re.I)
    return text


def final_gate(words: list[str]) -> int:
    """Scan everything that would be published; non-zero exit on any hit."""
    bad = []
    for f in REPO.rglob("*"):
        relative = f.relative_to(REPO)
        if not f.is_file():
            continue
        if any(part in JUNK_DIRS for part in relative.parts[:-1]):
            continue
        if f.suffix.lower() in JUNK_EXTS:
            continue
        if publication_hit(f.read_bytes(), words):
            bad.append(relative.as_posix())
    if bad:
        print("PUBLICATION GATE FAILED — do not publish; offending files:")
        for b in bad:
            print("   ", b)
        return 1
    print("publication gate: clean")
    return 0


# ---------------------------------------------------------------- progress
def newest_progress_block(
    project: dict,
    progress: Path | None = None,
) -> tuple[str, str]:
    """Return the newest matching block from the project's domain ledger.

    Precedence (each PROGRESS.md is append-order newest-first, but campaigns
    interleave, so pure keyword-first-match is unreliable):
      1. block containing the project's anchor regex (if set)
      2. block whose header line contains the slug (word-boundary)
      3. block whose body contains the slug as a path ("<slug>/")
      4. first block matching the keyword list
    """
    if progress is None:
        domain = DOMAINS[project_domain(project)]
        progress_name = domain["progress"]
        if not progress_name:
            return "", ""
        progress = domain["root"] / progress_name
    if not progress.exists():
        return "", ""
    lines = progress.read_text(errors="replace").splitlines()
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
            if len(text) >= 300 and re.search(re.escape(anchor), text, re.I):
                return hdr, text[:6500]
    for _, hdr, bd in blocks:
        text = "\n".join([hdr] + bd)
        if (len(text) >= 300
                and re.search(
                    rf"(?<![a-z0-9]){re.escape(slug)}(?![a-z0-9])",
                    hdr.lower(),
                )):
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


def quant_trading_summary(words: list[str]) -> tuple[str, str, str]:
    """Build an allowlisted public summary from the quant-trading registries."""
    root = DOMAINS["quant-trading"]["root"]
    workspace = json.loads((root / "workspace.json").read_text())
    registry = json.loads(
        (root / "explorations" / "strategy-registry.json").read_text()
    )
    summary = registry["summary"]
    purpose = redact(str(workspace["purpose"]), words)
    status = redact(str(workspace["status"]), words)
    conclusion = redact(str(summary["final_conclusion"]), words)
    updated = redact(str(registry["updated_at"]), words)
    markdown = f"""# Quantitative trading research

**Status (from `workspace.json`):** {status}

{purpose}

## Latest registered conclusion

{conclusion}

## Registered campaign counts

- Independent wave-one directions: {summary["independent_wave_one_directions"]}
- Wave-one selection trials: {summary["wave_one_selection_trials"]}
- Initial validation passes: {summary["initial_validation_passes"]}
- Final validated or deployable candidates: {summary["final_validated_or_deployable_candidates"]}
- Focused directions: {summary["focused_directions"]}

Registry updated: `{updated}`.

Only allowlisted status fields are rendered here. Raw market data, CSV return
series, active capture-process names, machine-local paths, and provisional
output trees are not published.
"""
    return purpose, conclusion[:180], markdown


# ---------------------------------------------------------------- mirror


def collect_sources(
    project: dict,
    source: Path | None = None,
) -> dict[str, Path]:
    """Return public relative paths mapped to source files."""
    src = source if source is not None else project_source(project)
    out: dict[str, Path] = {}
    for rel in project.get("papers", ()):
        f = src / rel
        if not f.is_file():
            continue
        if f.stat().st_size > LIMIT_KB * 1024:
            print(f"  skip (>{LIMIT_KB}KB): {project['slug']}/{rel}")
            continue
        out[Path(rel).as_posix()] = f
    if not project.get("walk", True) or not src.is_dir():
        return out

    exclude_dirs = set(project.get("exclude_dirs", ()))
    exclude_paths = tuple(
        Path(path).as_posix().strip("/")
        for path in project.get("exclude_paths", ())
    )
    include_exts = {
        suffix.lower() for suffix in project.get("include_exts", ())
    }
    for f in sorted(src.rglob("*")):
        if not f.is_file():
            continue
        rel = f.relative_to(src)
        rel_text = rel.as_posix()
        if any(p in JUNK_DIRS or p in exclude_dirs
               for p in rel.parts[:-1]):
            continue
        if any(rel_text == prefix or rel_text.startswith(prefix + "/")
               for prefix in exclude_paths):
            continue
        if f.suffix.lower() in JUNK_EXTS or f.name == ".DS_Store":
            continue
        if include_exts and f.suffix.lower() not in include_exts:
            continue
        if f.stat().st_size > LIMIT_KB * 1024:
            continue
        out.setdefault(rel_text, f)
    return out


def sync_materials(project: dict, words: list[str]) -> tuple[list[dict], dict]:
    """Incrementally mirror one project; return file metadata and statistics."""
    dest = REPO / project_material(project)
    dest.mkdir(parents=True, exist_ok=True)
    files: list[dict] = []
    stats = dict(written=0, unchanged=0, pruned=0, withheld=0)
    kept: set[str] = set()
    for rel, srcf in collect_sources(project).items():
        data = srcf.read_bytes()
        if source_blocked(project, data, words):
            stats["withheld"] += 1
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
    # Prune files that no longer belong to this project's public mirror.
    for f in sorted(dest.rglob("*")):
        if f.is_file() and f.relative_to(dest).as_posix() not in kept:
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
      var config = (window.__LINKS__ || {})[el.id];
      if (!config) return;
      var base = typeof config === 'string' ? config : config.base;
      var root = typeof config === 'string' ? null : config.root;
      var nodes = el.querySelectorAll('a[href], img[src]');
      for (var i = 0; i < nodes.length; i++) {
        var n = nodes[i];
        var attr = n.tagName === 'IMG' ? 'src' : 'href';
        var h = n.getAttribute(attr);
        if (!h) continue;
        if (/^(?:https?:|mailto:|tel:|#|\\/\\/|\\/)/i.test(h)) continue;
        var prefix = root && /^(?:\\.\\/)?(?:README|RESULTS|PROGRESS)\\.md(?:[#?].*)?$/i.test(h)
          ? root : base;
        n.setAttribute(attr, prefix + h);
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
         links: dict[str, object] | None = None) -> str:
    md_json = json.dumps(md_map, ensure_ascii=False).replace("<", "\\u003c")
    lk_json = json.dumps(links or {}, ensure_ascii=False)
    home = "../" * depth + "index.html"
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
<hr><p class="tag">Auto-generated by tools/sync.py · public research mirror</p>
</body></html>"""


def build_problem(
    e: dict,
) -> tuple[str, dict[str, str], dict[str, object]]:
    mat = f"../{e['material']}/"
    groups: dict[str, list[dict]] = {}
    for f in e["files"]:
        top = f["path"].split("/", 1)[0] if "/" in f["path"] else "(root)"
        groups.setdefault(top, []).append(f)
    parts = []
    for group in sorted(groups, key=lambda key: (key != "(root)", key.lower())):
        files = groups[group]
        mb = sum(item["size"] for item in files) / 1e6
        items = "".join(
            f'<li><a href="{mat}{esc(item["path"])}">{esc(item["path"])}</a>'
            f' <span class="fsize">{max(item["size"] // 1024, 1)} KB</span></li>'
            for item in files
        )
        label = "project root" if group == "(root)" else group + "/"
        opened = " open" if group == "(root)" else ""
        parts.append(
            f'<details{opened}><summary><code>{esc(label)}</code>'
            f' — {len(files)} files, {mb:.1f} MB</summary>'
            f'<ul class="files">{items}</ul></details>'
        )
    files_html = "\n".join(parts) or '<p class="tip">No artifacts copied.</p>'

    domain = e["domain"]
    domain_config = DOMAINS[domain]
    shared_links = " · ".join(
        f'<a href="../materials/domains/{esc(domain)}/{esc(name)}">'
        f"<code>{esc(name)}</code></a>"
        for name in e["shared_files"]
    )
    shared_html = (
        f'<p class="tag">Domain documents: {shared_links}</p>'
        if shared_links else ""
    )

    md_map: dict[str, str] = {}
    if e["progress"]:
        md_map["md-progress"] = e["progress"]
        prog = (
            '<h2>Latest progress</h2><p class="tag">newest matching entry, '
            f'quoted from <code>{esc(domain)}/PROGRESS.md</code></p>'
            '<div class="markdown-body" id="md-progress"></div>'
        )
    elif domain_config["progress"]:
        prog = '<p class="tip">No PROGRESS.md entry matched this target.</p>'
    else:
        prog = ""
    if e["readme_head"]:
        md_map["md-readme"] = e["readme_head"]
        if e.get("summary"):
            heading = "Research status"
            note = "allowlisted fields from the canonical JSON registries"
        else:
            heading = "Target README"
            note = "first lines, verbatim"
        readme = (
            f"<h2>{heading}</h2><p class=\"tag\">{note}</p>"
            '<div class="markdown-body" id="md-readme"></div>'
        )
    else:
        readme = ""
    session = (
        f' · session <code>{esc(e["session"])}</code>'
        if e["session"] else ""
    )
    body = f"""<h1>{esc(e['name'])}</h1>
<p class="tag">source <code>{esc(e['source'])}</code>{session} ·
last change {esc(e['updated'])}</p>
{shared_html}
{prog}
{readme}
<h2>Code &amp; artifacts</h2>
<p class="tag">mirrored from the source tree; oversized, raw, provisional,
third-party, and publication-gated files are excluded</p>
{files_html}"""
    progress_material = (
        "../materials/" if domain == "math"
        else f"../materials/{domain}/"
    )
    return body, md_map, {
        "md-progress": {
            "base": progress_material,
            "root": f"../materials/domains/{domain}/",
        },
        "md-readme": mat,
    }


def build_index(
    entries: list[dict],
    session_md: str,
    shared_files: dict[str, list[str]],
) -> tuple[str, dict[str, str], dict[str, str]]:
    sections = []
    for domain, config in DOMAINS.items():
        domain_entries = sorted(
            (entry for entry in entries if entry["domain"] == domain),
            key=lambda entry: entry["name"].lower(),
        )
        if not domain_entries:
            continue
        rows = []
        for entry in domain_entries:
            rows.append(
                "<tr>"
                f'<td><a href="{entry["page"]}">{esc(entry["name"])}</a></td>'
                f'<td><code>{esc(entry["source"])}</code></td>'
                f'<td>{esc(entry["updated"][:10])}</td>'
                f'<td>{esc(strip_md(entry["headline"]))}</td>'
                f'<td>{len(entry["files"])}</td>'
                "</tr>"
            )
        links = " · ".join(
            f'<a href="materials/domains/{esc(domain)}/{esc(name)}">'
            f"<code>{esc(name)}</code></a>"
            for name in shared_files.get(domain, [])
        )
        root_docs = f'<p class="tag">Root documents: {links}</p>' if links else ""
        sections.append(
            f"<h2>{esc(config['name'])}</h2>"
            f"<p class=\"tag\">source <code>~/jinleic-workspace/{esc(domain)}</code></p>"
            f"{root_docs}"
            "<table>"
            "<tr><th>Target</th><th>source</th><th>last change</th>"
            "<th>latest headline</th><th>files</th></tr>"
            f"{''.join(rows)}"
            "</table>"
        )
    body = f"""<h1>what-is-this — research progress tracker</h1>
<p class="tag">Public mirror of the tracked research under
<code>math/</code>, <code>physics/</code>, <code>cs/</code>, and
<code>quant-trading/</code>. Target pages quote their README and newest matching
domain progress entry, then link core code and small result artifacts. Raw data,
provisional trees, oversized files, and publication-gated content stay private.
Writes are incremental.</p>
<p class="tip">There is no separate <code>quantum/</code> source tree.
Quantum-computing targets live under <code>physics/</code>; algebraic quantum
LDPC work also lives at <code>math/qec/</code>.</p>

{''.join(sections)}

<h2>How to update</h2>
<pre>cd ~/jinleic-workspace/repos/what-is-this
python3 tools/sync.py        # incremental; fails closed on the publication gate
git add -A && git commit -m "sync $(date -u +%F)" && git push</pre>

<h2>Recorded mathematics sessions</h2>
<div class="markdown-body" id="md-sessions"></div>"""
    return body, {"md-sessions": session_md}, {"md-sessions": ""}


def build_readme(
    entries: list[dict],
    shared_files: dict[str, list[str]],
) -> str:
    sections = []
    for domain, config in DOMAINS.items():
        domain_entries = sorted(
            (entry for entry in entries if entry["domain"] == domain),
            key=lambda entry: entry["name"].lower(),
        )
        if not domain_entries:
            continue
        rows = []
        for entry in domain_entries:
            headline = strip_md(entry["headline"]).replace("|", "\\|")
            rows.append(
                f"| [{entry['name']}]({entry['page']}) "
                f"| `{entry['source']}` | {entry['updated'][:10]} | {headline} |"
            )
        root_links = ", ".join(
            f"[`{name}`](materials/domains/{domain}/{name})"
            for name in shared_files.get(domain, [])
        )
        root_line = f"\nRoot documents: {root_links}.\n" if root_links else ""
        sections.append(
            f"## {config['name']}\n\n"
            f"Source: `~/jinleic-workspace/{domain}`.\n"
            f"{root_line}\n"
            "| Target | source | last change | latest headline |\n"
            "|---|---|---|---|\n"
            + "\n".join(rows)
        )
    sections_text = "\n\n".join(sections)
    return f"""# what-is-this — research progress tracker

Public progress tracker for research under `math/`, `physics/`, `cs/`, and
`quant-trading/`. Quantum-computing targets are under `physics/`; algebraic
quantum LDPC work is also under `math/qec/`. There is no separate top-level
`quantum/` source directory.

Regenerated by `tools/sync.py`: target README and progress excerpts, domain root
results/progress documents, core code, and small artifacts. Raw data,
provisional or third-party trees, oversized files, private keywords, and
high-confidence secret shapes are excluded. The quantitative-trading mirror
also withholds machine-local paths. Writes are incremental.

{sections_text}

## Update

```sh
python3 tools/sync.py
git add -A && git commit -m "sync $(date -u +%F)"
git push
```

`session-map.md` records the existing mathematics session mapping.
`tools/sync.py` is the source of truth for mirror paths, size caps, exclusions,
and the fail-closed publication gate.
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
| `01a03c0c-c36d-73b7-96b7-f1203909eb76` | 2026-08-25 | Erdős Problem #389 (e389) |
| `01a03b5f-2600-7019-adf6-b03bd015c2ce` | 2026-08-25 | UC audit + Liu Hypothesis 1 submission package (liu_h1) |

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
    prev_domains = prev.get("domains", {})

    total_bytes = 0
    pages_rewritten = 0
    shared_files: dict[str, list[str]] = {}
    domain_manifest: dict[str, dict] = {}
    for domain, config in DOMAINS.items():
        names = config["shared_files"]
        if not names:
            shared_files[domain] = []
            continue
        mirror = dict(
            domain=domain,
            slug="_domain",
            dirname=".",
            papers=names,
            walk=False,
            material=f"materials/domains/{domain}",
        )
        files, stats = sync_materials(mirror, words)
        shared_files[domain] = [item["path"] for item in files]
        byte_count = sum(item["size"] for item in files)
        total_bytes += byte_count
        sig = hashlib.sha256(
            json.dumps(files, sort_keys=True).encode()
        ).hexdigest()[:16]
        previous = prev_domains.get(domain, {})
        updated = (
            previous.get("updated", ts)
            if previous.get("sig") == sig else ts
        )
        domain_manifest[domain] = dict(
            updated=updated,
            sig=sig,
            bytes=byte_count,
            files=files,
        )
        print(
            f"{domain + '/_domain':24s} files={len(files)} "
            f"({byte_count / 1e6:.1f} MB): "
            f"{stats['written']} written, {stats['unchanged']} unchanged, "
            f"{stats['pruned']} pruned, {stats['withheld']} withheld"
        )

    entries: list[dict] = []
    for project in PROJECTS:
        domain = project_domain(project)
        source = project_source(project)
        if project.get("summary") == "quant-trading":
            tagline, headline, readme_head = quant_trading_summary(words)
            header = block = ""
        else:
            readme_path = source / "README.md"
            readme_head = redact(read_head(readme_path), words)
            tagline = redact(readme_tagline(readme_path), words)
            header, block = newest_progress_block(project)
            header = redact(header, words)
            block = redact(block, words)
            match = re.search(r"(?m)^#{1,3}\s+.*$", block)
            headline = (
                match.group(0) if match else first_line_nonempty(block)
            )[:180]
            if not headline:
                headline = strip_md(tagline)[:180]

        files, stats = sync_materials(project, words)
        byte_count = sum(item["size"] for item in files)
        total_bytes += byte_count
        sig = hashlib.sha256(
            (block + "\x00" + readme_head + "\x00"
             + json.dumps(files, sort_keys=True)).encode()
        ).hexdigest()[:16]
        key = project_manifest_key(project)
        previous = prev_projects.get(key, {})
        updated = (
            previous.get("updated", ts)
            if previous.get("sig") == sig else ts
        )
        entry = dict(
            id=key,
            slug=project["slug"],
            domain=domain,
            name=project["name"],
            dirname=project["dirname"],
            source=project_source_label(project),
            material=project_material(project),
            session=project.get("session", ""),
            summary=project.get("summary"),
            tagline=tagline,
            headline=headline,
            progress_header=header,
            progress=block,
            readme_head=readme_head,
            files=files,
            shared_files=shared_files.get(domain, []),
            updated=updated,
            sig=sig,
            page=project_page(project),
        )
        entries.append(entry)
        body, md_map, links = build_problem(entry)
        page_path = REPO / entry["page"]
        depth = len(page_path.relative_to(REPO).parent.parts)
        if write_if_changed(
            page_path,
            page(project["name"], body, md_map, depth=depth, links=links),
        ):
            pages_rewritten += 1
        print(f"{key:24s} {headline[:100] or '(no entry)'}")
        print(
            f"{'':24s} files={len(files)} ({byte_count / 1e6:.1f} MB): "
            f"{stats['written']} written, {stats['unchanged']} unchanged, "
            f"{stats['pruned']} pruned, {stats['withheld']} withheld"
        )

    domain_rank = {name: rank for rank, name in enumerate(DOMAINS)}
    entries_sorted = sorted(
        entries,
        key=lambda entry: (
            domain_rank[entry["domain"]],
            entry["name"].lower(),
        ),
    )
    session_md = build_session_map()
    write_if_changed(REPO / "session-map.md", session_md)
    write_if_changed(
        REPO / "README.md",
        build_readme(entries_sorted, shared_files),
    )
    idx_body, idx_md, idx_links = build_index(
        entries_sorted,
        session_md,
        shared_files,
    )
    if write_if_changed(
        REPO / "index.html",
        page(
            "what-is-this — research progress tracker",
            idx_body,
            idx_md,
            depth=0,
            links=idx_links,
        ),
    ):
        pages_rewritten += 1

    update_times = [entry["updated"] for entry in entries]
    update_times.extend(
        domain["updated"] for domain in domain_manifest.values()
    )
    manifest = dict(
        updated=max(update_times, default=ts),
        limit_kb=LIMIT_KB,
        domains=domain_manifest,
        projects={
            entry["id"]: dict(
                domain=entry["domain"],
                updated=entry["updated"],
                sig=entry["sig"],
                page=entry["page"],
                material=entry["material"],
                headline=entry["headline"],
                progress_header=entry["progress_header"],
                bytes=sum(item["size"] for item in entry["files"]),
                files=entry["files"],
            )
            for entry in entries
        },
    )
    write_if_changed(
        sj,
        json.dumps(manifest, indent=1, sort_keys=True, ensure_ascii=False),
    )

    print(
        f"\ntotal mirrored: {total_bytes / 1e6:.1f} MB; "
        f"pages rewritten this run: {pages_rewritten}"
    )
    rc = final_gate(words)
    if rc == 0:
        print("clean — sync with: git add -A && git commit && git push")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
