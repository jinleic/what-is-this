#!/usr/bin/env python3
"""Stage the h10q deliverable bundle and zip it.

Layout (dated deliverable tree):
    deliverables/<DATE>/papers/      main + companion drafts, sections/
    deliverables/<DATE>/artifacts/   data artifacts cited by the papers
    deliverables/<DATE>/scripts/     checked-in regenerators / replay scripts
    deliverables/<DATE>/ledger/      THEOREMS / RESULTS / NOTES / README / CONDITIONAL
    deliverables/<DATE>/MANIFEST.md  human-readable inventory + provenance labels
    deliverables/<DATE>/SHA256SUMS   checksums of every staged file

The ZIP is written OUTSIDE the staged tree (deliverables/h10q-bundle-<DATE>.zip)
so it can never contain itself.

Run:  python3 deliverables/make_bundle.py           (from math/h10q)
"""
from __future__ import annotations

import hashlib
import json
import os
import zipfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # math/h10q
DELIV = ROOT / "deliverables"
DATE = os.environ.get("H10Q_BUNDLE_DATE", "2026-08-19")
STAGE = DELIV / DATE
ZIP_PATH = DELIV / f"h10q-bundle-{DATE}.zip"

# (source, staged subdir, provenance label)
REGENERABLE = "repo-regenerable"
EVIDENCE = "persisted evidence (no checked-in producer)"
AUTHORITY = "frozen authority (suite-asserted)"

ARTIFACTS = [
    ("data/l17_matched.jsonl", REGENERABLE),
    ("data/l17_badroots_closures.jsonl", REGENERABLE),
    ("data/l17_classfactors.jsonl", REGENERABLE),
    ("data/l17_schinzel_audit.jsonl", REGENERABLE),
    ("data/l17_horizon101.jsonl", REGENERABLE),
    ("data/l17_horizon103.jsonl", REGENERABLE),
    ("data/l17_horizon107.jsonl", REGENERABLE),
    ("data/l17_horizon109.jsonl", REGENERABLE),
    ("data/l17_horizon109_crossk.jsonl", REGENERABLE),
    ("data/l17_horizon113.jsonl", REGENERABLE),
    ("data/l17_horizon127.jsonl", REGENERABLE),
    ("data/l17_horizon131.jsonl", REGENERABLE),
    ("data/l17_horizon137.jsonl", REGENERABLE),
    ("data/l17_composite_w.jsonl", REGENERABLE),
    ("data/l17_stepii.jsonl", REGENERABLE),
    ("data/l18_schinzel_implies_h.jsonl", REGENERABLE),
    ("data/l18_fivewall.jsonl", REGENERABLE),
    ("data/l18_route131.jsonl", REGENERABLE),
    ("data/l18_horizon_sweep_a.jsonl", REGENERABLE),
    ("data/l18_horizon_sweep_b.jsonl", REGENERABLE),
    ("data/l18_divergence_model.jsonl", REGENERABLE),
    ("data/l13h_all_closures.json", AUTHORITY),
    ("data/l6_witnesses.jsonl", AUTHORITY),
    ("data/l9_steered.jsonl", AUTHORITY),
    ("data/l15_remainders.jsonl", EVIDENCE),
    ("data/l15_density.json", EVIDENCE),
    ("data/l16_char.jsonl", EVIDENCE),
    ("data/l17_sieve.jsonl", EVIDENCE),
    ("data/l17_badroots.jsonl", EVIDENCE),
    ("data/l17_ratemodel_censored.jsonl", EVIDENCE),
    ("data/l17_cofactor.jsonl", EVIDENCE),
    ("data/l17_cofactor_nonclosure.jsonl", EVIDENCE),
    ("data/litscout_h10q.md", EVIDENCE),
]

SCRIPTS = [
    "h10q.py",
    "l13_filter.py",
    "l12_class.py",
    "l13h_scan.py",
    "l17_ratemodel_padic.py",
    "l17_badroots_closures.py",
    "l17_matched.py",
    "l17_classrisk.py",
    "l17_schinzel_audit.py",
    "l17_horizon101.py",
    "l17_horizon103.py",
    "l17_horizon107.py",
    "l17_horizon109.py",
    "l17_horizon109_crossk.py",
    "l17_horizon113.py",
    "l17_horizon127.py",
    "l17_horizon131.py",
    "l17_horizon137.py",
    "l17_composite_w.py",
    "l18_schinzel_implies_h.py",
    "l18_fivewall.py",
    "l18_route131.py",
    "l18_horizon_sweep_a.py",
    "l18_horizon_sweep_b.py",
    "l18_divergence_model.py",
    "l16_emergent.py",
    "l14_replay_all.py",
]

LEDGER = ["THEOREMS.md", "RESULTS.md", "NOTES.md", "README.md", "CONDITIONAL.md"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stage_file(src: Path, dest: Path) -> bool:
    if not src.is_file():
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(src.read_bytes())
    return True


def main() -> None:
    staged: list[tuple[str, str, int, str]] = []   # relpath, label, size, sha
    missing: list[str] = []

    for rel, label in ARTIFACTS:
        src = ROOT / rel
        dest = STAGE / "artifacts" / Path(rel).name
        if stage_file(src, dest):
            staged.append((f"artifacts/{Path(rel).name}", label, dest.stat().st_size, sha256(dest)))
        else:
            missing.append(rel)

    for rel in SCRIPTS:
        src = ROOT / rel
        dest = STAGE / "scripts" / rel
        if stage_file(src, dest):
            staged.append((f"scripts/{rel}", "generator/replay script", dest.stat().st_size, sha256(dest)))
        else:
            missing.append(rel)

    for rel in LEDGER:
        src = ROOT / rel
        dest = STAGE / "ledger" / rel
        if stage_file(src, dest):
            staged.append((f"ledger/{rel}", "project ledger", dest.stat().st_size, sha256(dest)))
        else:
            missing.append(rel)

    for tex in sorted((STAGE / "papers").rglob("*.tex")):
        rel = tex.relative_to(STAGE).as_posix()
        staged.append((rel, "paper draft", tex.stat().st_size, sha256(tex)))

    # SHA256SUMS
    sums = STAGE / "SHA256SUMS"
    sums.write_text("".join(f"{sha}  {rel}\n" for rel, _, _, sha in sorted(staged)))

    # MANIFEST
    lines = [
        f"# h10q deliverable bundle --- {DATE}",
        "",
        "Two paper drafts plus the artifacts, generators and ledger they cite.",
        "",
        "## Provenance labels",
        "",
        "- **repo-regenerable**: a checked-in script in `scripts/` regenerates the artifact.",
        "- **frozen authority (suite-asserted)**: byte-identical match is enforced by the test suites.",
        "- **persisted evidence (no checked-in producer)**: reviewed data whose original",
        "  session-side generator was not recovered into the repository.",
        "- **generator/replay script**, **paper draft**, **project ledger**: as named.",
        "",
        "## Contents",
        "",
        "| file | provenance | bytes |",
        "|---|---|---|",
    ]
    for rel, label, size, _ in sorted(staged):
        lines.append(f"| `{rel}` | {label} | {size} |")
    lines += [
        "",
        "## Verification",
        "",
        "```sh",
        "shasum -a 256 -c SHA256SUMS      # from inside this directory",
        "```",
        "",
        "## Reproduction (repo-regenerable chain, run from math/h10q)",
        "",
        "```sh",
        "python3 l17_ratemodel_padic.py            # p-adic rate model / bad-mass engine",
        "python3 l17_badroots_closures.py          # closure-authority bad-root tables",
        "python3 l17_matched.py                    # matched statistics artifact",
        "python3 l17_classrisk.py                  # per-class window factors",
        "python3 l17_schinzel_audit.py             # Schinzel local-condition audit",
        "python3 l17_horizon101.py                 # off-grid closure replays (103/107/109 likewise)",
        "python3 l17_horizon109_crossk.py          # ~17 min deterministic cross-k scan",
        "python3 h10q.py                           # default suite (exit 0 required)",
        "python3 h10q.py --extended                # extended suite (exit 0 required)",
        "```",
        "",
    ]
    if missing:
        lines += ["## Missing at build time", ""] + [f"- `{m}`" for m in missing] + [""]
    (STAGE / "MANIFEST.md").write_text("\n".join(lines))

    # ZIP (outside the staged tree)
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(STAGE.rglob("*")):
            if path.is_file():
                zf.write(path, arcname=f"h10q-{DATE}/{path.relative_to(STAGE).as_posix()}")

    print(json.dumps({
        "staged_files": len(staged) + 2,
        "missing": missing,
        "zip": str(ZIP_PATH.relative_to(ROOT)),
        "zip_bytes": ZIP_PATH.stat().st_size,
        "zip_sha256": sha256(ZIP_PATH),
    }, indent=1))


if __name__ == "__main__":
    main()
