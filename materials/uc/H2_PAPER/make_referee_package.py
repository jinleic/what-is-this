#!/usr/bin/env python3
"""Build the external referee package for c' = 1 - m* and c'' = 1 - m_{16/25}.

The package is DERIVED, never stored in the repository: this script is the single
authority for what a referee needs, and it copies the canonical files from their
one home under math/ into a distribution directory plus a deterministic tarball.
Nothing here is a second copy that can drift -- rerun the script and any change
in a canonical file shows up as a changed hash in the generated MANIFEST.

The file list below is closed under three relations, each verified by
`--check` before anything is written:

  * everything uc/H2_PAPER/main.tex names via \\path{...};
  * the eight CHAIN links of uc/verification/replay_h2_chain.py (module +
    artifact) and the 24 hashes they pin;
  * the transitive local imports of those eight modules, which reach exactly
    three further modules (liu9_objective.py, liu9_binding.py,
    liu9_cprime_frontier.py) and one further artifact
    (liu9-cprime-frontier.json, asserted at run time by the eighth module).

PROGRESS.md is included because the driver's check [A] asserts all 24 pinned
hashes against it: shipping the ledger is what makes the anti-drift check
runnable, not decoration.

Standard library only.  Run from math/:

    ./.venv/bin/python -I -B uc/H2_PAPER/make_referee_package.py
    ./.venv/bin/python -I -B uc/H2_PAPER/make_referee_package.py --check
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import io
import json
import re
import shutil
import sys
import tarfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

HERE = Path(__file__).resolve().parent          # uc/H2_PAPER
ROOT = HERE.parent.parent                       # math/
DIST = HERE / "dist"

MANUSCRIPT = "uc/H2_PAPER/main.tex"
DRIVER = "uc/verification/replay_h2_chain.py"
LEDGER = "PROGRESS.md"

# Pinned helper modules: imported by the chain but not themselves chain links.
HELPER_MODULES = (
    "uc/liu9_objective.py",
    "uc/liu9_binding.py",
    "uc/liu9_cprime_frontier.py",
)
# Pinned input artifact of the eighth link, asserted by it at run time.
HELPER_ARTIFACTS = ("uc/verification/results/liu9-cprime-frontier.json",)

EXTRAS = (
    ("requirements-freeze.txt", "third-party version lock; byte-identical replay depends on it"),
    ("uc/H2_PAPER/REFEREE.md", "referee entry point: replay recipe and independent-verification checklist"),
    ("uc/H2_PAPER/make_referee_package.py", "this generator; the package's own provenance"),
    ("uc/LITERATURE_ORIGINALITY.md",
     "originality audit; single authority for the literature position of c' and c''"),
    ("uc/literature/search_protocol.json",
     "the queried endpoints and their coverage limits"),
    ("uc/literature/arxiv_refresh_2026-09-03.json",
     "arXiv API record behind the 2026-09-03 literature claim"),
    ("LIU_H1/literature/pdfs/liu_2023_arxiv_2306.08824v1.pdf",
     "Liu 2023, the cited primary source; required only to re-run uc/liu9_cprime_frontier.py itself"),
)

ROLES = {
    MANUSCRIPT: "manuscript: both constants, every machine-assisted theorem citing its artifact",
    LEDGER: "ledger; driver check [A] asserts all 24 pinned hashes against it",
    DRIVER: "replay driver: 8 links x 2 runs, 155 checks",
    "uc/liu9_objective.py": "pinned helper: audited transcription of Liu's objective",
    "uc/liu9_binding.py": "pinned helper: solver for Liu's equations (87)-(90)",
    "uc/liu9_cprime_frontier.py": "pinned helper: the exact candidate balls for the scaled protocol",
    "uc/verification/results/liu9-cprime-frontier.json": "pinned input of link 8, asserted at its run time",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_chain() -> List[dict]:
    """Read the CHAIN table out of the driver by AST, without importing it."""
    source = (ROOT / DRIVER).read_text(encoding="utf-8")
    tree = ast.parse(source)
    fields = None
    chain = None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == "Link":
            fields = [s.target.id for s in node.body if isinstance(s, ast.AnnAssign)]
        if isinstance(node, ast.AnnAssign) and getattr(node.target, "id", None) == "CHAIN":
            chain = node.value
    if fields is None or chain is None:
        raise SystemExit("could not locate Link/CHAIN in the driver")
    links = []
    for call in chain.elts:
        values = [ast.literal_eval(a) for a in call.args]
        links.append(dict(zip(fields, values)))
    return links


def local_imports(path: Path) -> List[str]:
    """Local (repo) modules imported by `path`, by AST; lazy imports included."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names: List[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names += [a.name.split(".")[0] for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names.append(node.module.split(".")[0])
    out = []
    for n in sorted(set(names)):
        candidate = ROOT / "uc" / (n + ".py")
        if candidate.is_file():
            out.append("uc/%s.py" % n)
    return out


def import_closure(seeds: Iterable[str]) -> List[str]:
    seen: set[str] = set()
    stack = list(seeds)
    while stack:
        rel = stack.pop()
        if rel in seen:
            continue
        seen.add(rel)
        stack += local_imports(ROOT / rel)
    return sorted(seen)


def manuscript_paths() -> List[str]:
    """Every repo path main.tex names via \\path{...} or \\hash{...}."""
    tex = (ROOT / MANUSCRIPT).read_text(encoding="utf-8")
    raw = re.findall(r"\\path\{([^}]*)\}", tex)
    out = set()
    for token in raw:
        token = token.replace("\\_", "_").strip()
        if not token or token.endswith("/") or "." not in token:
            continue
        for candidate in (token, "uc/" + token, "uc/verification/results/" + token):
            if (ROOT / candidate).is_file():
                out.add(candidate)
                break
    return sorted(out)


def build_file_list() -> Tuple[List[Tuple[str, str]], List[dict]]:
    links = load_chain()
    modules = [l["module"] for l in links]
    artifacts = [l["artifact"] for l in links]
    ordered: List[Tuple[str, str]] = [
        (MANUSCRIPT, ROLES[MANUSCRIPT]),
        (LEDGER, ROLES[LEDGER]),
        (DRIVER, ROLES[DRIVER]),
    ]
    for i, l in enumerate(links, 1):
        ordered.append((l["module"], "chain module %d (%s)" % (i, l["name"])))
    for rel in HELPER_MODULES:
        ordered.append((rel, ROLES[rel]))
    for i, l in enumerate(links, 1):
        ordered.append((l["artifact"], "artifact %d (%s), claim_status %s, digest scope %s"
                        % (i, l["name"], l["claim_status"], l["digest_scope"])))
    for rel in HELPER_ARTIFACTS:
        ordered.append((rel, ROLES[rel]))
    for rel, role in EXTRAS:
        ordered.append((rel, role))
    seen = set()
    unique = []
    for rel, role in ordered:
        if rel not in seen:
            seen.add(rel)
            unique.append((rel, role))
    return unique, links


def check(files: List[Tuple[str, str]], links: List[dict]) -> List[str]:
    """Closure and pin checks.  Returns a list of failures."""
    listed = {rel for rel, _ in files}
    failures: List[str] = []

    for rel, _ in files:
        if not (ROOT / rel).is_file():
            failures.append("missing file: %s" % rel)

    # 1. the 24 chain pins match disk and occur in the ledger
    ledger = (ROOT / LEDGER).read_text(encoding="utf-8")
    pins = set()
    for l in links:
        pins |= {l["file_sha256"], l["report_sha256"], l["tool_sha256"]}
        if sha256(ROOT / l["artifact"]) != l["file_sha256"]:
            failures.append("artifact hash drift: %s" % l["artifact"])
        if sha256(ROOT / l["module"]) != l["tool_sha256"]:
            failures.append("module hash drift: %s" % l["module"])
        body = json.loads((ROOT / l["artifact"]).read_text(encoding="utf-8"))
        stated = body.get("report_sha256")
        if l["digest_scope"] == "empty":
            body["report_sha256"] = ""
        else:
            body.pop("report_sha256", None)
        blob = json.dumps(body, sort_keys=True, separators=(",", ":"))
        if l["digest_scope"] == "omit+nl":
            blob += "\n"
        computed = hashlib.sha256(blob.encode("utf-8")).hexdigest()
        if computed != stated or stated != l["report_sha256"]:
            failures.append("internal digest drift: %s" % l["artifact"])
    if len(pins) != 24:
        failures.append("expected 24 distinct chain pins, found %d" % len(pins))
    for pin in sorted(pins):
        if pin not in ledger:
            failures.append("pin absent from %s: %s" % (LEDGER, pin[:16]))

    # 2. import closure adds nothing outside the list
    closure = import_closure([l["module"] for l in links])
    for rel in closure:
        if rel not in listed:
            failures.append("import closure needs an unlisted module: %s" % rel)
    # 3. every repo path the manuscript names is in the list
    for rel in manuscript_paths():
        if rel not in listed:
            failures.append("manuscript names an unlisted file: %s" % rel)

    return failures


def write_package(files: List[Tuple[str, str]], links: List[dict]) -> Tuple[Path, Path]:
    if DIST.exists():
        shutil.rmtree(DIST)
    payload = DIST / "referee-package"
    rows = []
    for rel, role in files:
        src = ROOT / rel
        dst = payload / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(src.read_bytes())
        rows.append((rel, sha256(src), src.stat().st_size, role))

    manifest = ["# MANIFEST -- generated by uc/H2_PAPER/make_referee_package.py; do not edit",
                "# path\tsha256\tbytes\trole"]
    manifest += ["%s\t%s\t%d\t%s" % r for r in rows]
    (payload / "MANIFEST.tsv").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    (payload / "SHA256SUMS").write_text(
        "".join("%s  %s\n" % (h, rel) for rel, h, _, _ in rows), encoding="utf-8")

    # deterministic tarball: sorted members, zeroed metadata, gzip mtime 0
    members = sorted([rel for rel, _, _, _ in rows] + ["MANIFEST.tsv", "SHA256SUMS"])
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w") as tar:
        for rel in members:
            path = payload / rel
            info = tarfile.TarInfo("referee-package/" + rel)
            data = path.read_bytes()
            info.size = len(data)
            info.mtime = 0
            info.mode = 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tar.addfile(info, io.BytesIO(data))
    import gzip
    tarball = DIST / "referee-package.tar.gz"
    with open(tarball, "wb") as fh:
        fh.write(gzip.compress(raw.getvalue(), compresslevel=9, mtime=0))
    return payload, tarball


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="verify closure and pins; write nothing")
    args = parser.parse_args(argv)

    files, links = build_file_list()
    failures = check(files, links)
    print("REFEREE_PACKAGE files=%d links=%d" % (len(files), len(links)))
    for f in failures:
        print("FAIL", f)
    if failures:
        print("REFEREE_PACKAGE FAIL (%d)" % len(failures))
        return 1
    print("closure: manuscript paths, 8 chain links, import closure, 24 ledger pins all OK")
    if args.check:
        print("REFEREE_PACKAGE CHECK PASS")
        return 0
    payload, tarball = write_package(files, links)
    total = sum((ROOT / rel).stat().st_size for rel, _ in files)
    print("payload  %s (%d bytes across %d files)" % (payload.relative_to(ROOT), total, len(files)))
    print("tarball  %s  sha256 %s" % (tarball.relative_to(ROOT), sha256(tarball)))
    print("REFEREE_PACKAGE PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
