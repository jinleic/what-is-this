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
    ("data/l19_tauzero.jsonl", REGENERABLE),
    ("data/l19_classexist.jsonl", REGENERABLE),
    ("data/l19_cell179.jsonl", REGENERABLE),
    ("data/l20_admissible.jsonl", REGENERABLE),
    ("data/l21_reducible_locus.jsonl", REGENERABLE),
    ("data/l21_irred_generic.jsonl", REGENERABLE),
    ("data/l21_irred_direct.jsonl", REGENERABLE),
    ("data/l22_elimination.jsonl", REGENERABLE),
    ("data/l22_reciprocal_cube.jsonl", REGENERABLE),
    ("data/l22_infinity.jsonl", REGENERABLE),
    ("data/l22_factor_tuple.jsonl", REGENERABLE),
    ("data/l22_square_branch.jsonl", REGENERABLE),
    ("data/l22_fiber_geometry.jsonl", REGENERABLE),
    ("data/l23_half_sieve.jsonl", REGENERABLE),
    ("data/l23_fibration.jsonl", REGENERABLE),
    ("data/l23_norm_section.jsonl", REGENERABLE),
    ("data/l23_rational_section.jsonl", REGENERABLE),
    ("data/l23_squareclass.jsonl", REGENERABLE),
    ("data/l23_multivar.jsonl", REGENERABLE),
    ("data/l23_absorption.jsonl", REGENERABLE),
    ("data/l24_diagonal_geometry.jsonl", REGENERABLE),
    ("data/l24_diagonal_arithmetic.jsonl", REGENERABLE),
    ("data/l24_diagonal_local.jsonl", REGENERABLE),
    ("data/l24_diagonal_search.jsonl", REGENERABLE),
    ("data/l25_scaled_coupling.jsonl", REGENERABLE),
    ("data/l26_reciprocal_tie.jsonl", REGENERABLE),
    ("data/l27_triangular_shear.jsonl", REGENERABLE),
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
    "l19_tauzero.py",
    "l19_classexist.py",
    "l19_cell179.py",
    "l20_admissible.py",
    "l21_reducible_locus.py",
    "l21_irred_generic.py",
    "l21_irred_direct.py",
    "l22_elimination.py",
    "l22_reciprocal_cube.py",
    "l22_infinity.py",
    "l22_factor_tuple.py",
    "l22_square_branch.py",
    "l22_fiber_geometry.py",
    "l23_half_sieve.py",
    "l23_fibration.py",
    "l23_norm_section.py",
    "l23_rational_section.py",
    "l23_squareclass.py",
    "l23_multivar.py",
    "l23_absorption.py",
    "l24_diagonal_geometry.py",
    "l24_diagonal_arithmetic.py",
    "l24_diagonal_local.py",
    "l24_diagonal_search.py",
    "l25_scaled_coupling.py",
    "l26_reciprocal_tie.py",
    "l27_triangular_shear.py",
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
        "python3 l19_tauzero.py                   # tau=0 wall criterion / class-side break",
        "python3 l19_classexist.py                # uniform aligned-class construction",
        "python3 l19_cell179.py                   # exact w=179 member closure",
        "python3 l20_admissible.py                # uniform admissibility audit",
        "nice -n 19 python3 l21_reducible_locus.py # explicit reducible locus / branch escape",
        "nice -n 19 python3 l21_irred_generic.py # generic and per-fiber irreducibility",
        "nice -n 19 python3 l21_irred_direct.py  # local no-go laws / reciprocal subfamily",
        "nice -n 19 python3 l22_elimination.py # PROVED fixed-Z theorem; scans EVIDENCE",
        "nice -n 19 python3 l22_reciprocal_cube.py # PROVED independent theorem; scans EVIDENCE",
        "nice -n 19 python3 l22_infinity.py # PROVED local structure; route OPEN",
        "nice -n 19 python3 l22_factor_tuple.py # PROVED criterion/no-go; compatibility OPEN",
        "nice -n 19 python3 l22_square_branch.py # PROVED free-lambda theorem; six-count OPEN/excluded",
        "nice -n 19 python3 l22_fiber_geometry.py # PROVED reductions; uniform lemma OPEN",
        "nice -n 19 python3 l23_half_sieve.py # PROVED small-prime-clean lower bound; full members OPEN",
        "nice -n 19 python3 l23_fibration.py # PROVED rank-10/Br and Capell barrier; member route OPEN",
        "nice -n 19 python3 l23_norm_section.py # PROVED scoped norm no-gos/reductions; global route OPEN",
        "nice -n 19 python3 l23_rational_section.py # PROVED no-section/genus-4 reductions; points OPEN",
        "nice -n 19 python3 l23_squareclass.py # PROVED dyadic/squareclass no-gos; scans EVIDENCE",
        "nice -n 19 python3 l23_multivar.py # PROVED tested-family degree bounds; no global closure",
        "nice -n 19 python3 l23_absorption.py # PROVED rank reductions; diagonal image PROVED empty",
        "nice -n 19 python3 l24_diagonal_geometry.py # PROVED both orientations Phi-empty; route FALSE",
        "nice -n 19 python3 l24_diagonal_arithmetic.py # PROVED exact dyadic closure; route CLOSED",
        "nice -n 19 python3 l24_diagonal_local.py # PROVED odd local theorem and dyadic supersession",
        "nice -n 19 python3 l24_diagonal_search.py # 4,026,282 zero hits EVIDENCE only",
        "python3 h10q.py                           # default suite (exit 0 required)",
        "python3 h10q.py --extended                # extended suite (exit 0 required)",
        "```",
        "",
    ]
    if missing:
        lines += ["## Missing at build time", ""] + [f"- `{m}`" for m in missing] + [""]
    (STAGE / "MANIFEST.md").write_text("\n".join(lines))

    # ZIP (outside the staged tree).  Write every member with the bundle date
    # rather than the source mtime: SHA256SUMS and MANIFEST.md are regenerated
    # on each run, so ZipFile.write() would otherwise make identical contents
    # produce a different archive hash.
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    year, month, day = (int(part) for part in DATE.split("-"))
    zip_timestamp = (year, month, day, 0, 0, 0)
    with zipfile.ZipFile(
        ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as zf:
        for path in sorted(STAGE.rglob("*")):
            if not path.is_file():
                continue
            arcname = f"h10q-{DATE}/{path.relative_to(STAGE).as_posix()}"
            info = zipfile.ZipInfo(arcname, date_time=zip_timestamp)
            info.create_system = 3
            info.external_attr = (path.stat().st_mode & 0xFFFF) << 16
            zf.writestr(
                info,
                path.read_bytes(),
                compress_type=zipfile.ZIP_DEFLATED,
                compresslevel=9,
            )

    print(json.dumps({
        "staged_files": len(staged) + 2,
        "missing": missing,
        "zip": str(ZIP_PATH.relative_to(ROOT)),
        "zip_bytes": ZIP_PATH.stat().st_size,
        "zip_sha256": sha256(ZIP_PATH),
    }, indent=1))


if __name__ == "__main__":
    main()
