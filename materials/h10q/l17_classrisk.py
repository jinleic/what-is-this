#!/usr/bin/env python3
"""Recompute the L17 p<=10000 class factors and closure-effort cross-check.

The output is deliberately a window statistic. It never uses
``small_factor_padic`` and never treats the residual singular-root bound as a
bound for primes above 10000.
"""
from __future__ import annotations

import json
import math
import os
import statistics
import time
from fractions import Fraction
from pathlib import Path

CAP = 10000
OUT_NAME = "l17_classfactors.jsonl"
REPORT = Path("/tmp/l17_classrisk.report")


def locate_root() -> Path:
    env = os.environ.get("H10Q_ROOT")
    candidates = []
    if env:
        candidates.append(Path(env))
    here = Path.cwd()
    candidates.extend((here, here / "math" / "h10q"))
    candidates.append(Path("/Users/jinleic/jinleic-workspace/math/h10q"))
    for candidate in candidates:
        if (candidate / "data" / "l17_badroots_closures.jsonl").exists():
            return candidate.resolve()
    raise FileNotFoundError("cannot locate math/h10q; set H10Q_ROOT")


ROOT = locate_root()
DATA = ROOT / "data"
BC_PATH = DATA / "l17_badroots_closures.jsonl"
CLOSURES_PATH = DATA / "l13h_all_closures.json"
OUT_PATH = DATA / OUT_NAME


def read_jsonl(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open() as handle:
        for i, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
            if i % 50 == 0:
                time.sleep(0.1)
    return rows


def as_fraction(value: object) -> Fraction:
    if isinstance(value, bool):
        return Fraction(int(value))
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, float):
        return Fraction(str(value))
    text = str(value)
    if "/" in text:
        num, den = text.split("/", 1)
        return Fraction(int(num), int(den))
    return Fraction(text)


def cell_key(cell: list) -> tuple[int, tuple[int, int]]:
    return int(cell[0]), (int(cell[1][0]), int(cell[1][1]))


def root_masses(class_row: dict) -> tuple[dict[int, Fraction], dict[int, Fraction], set[int]]:
    """Return known bad mass, unresolved residual upper mass, excluded p set."""
    known: dict[int, Fraction] = {}
    residual: dict[int, Fraction] = {}
    excluded: set[int] = set()
    for entry in class_row["bad_root_lists"]:
        p = int(entry["p"])
        if p > CAP:
            continue
        if entry.get("excluded"):
            excluded.add(p)
        for root in entry.get("roots", []):
            status = root.get("status")
            resolved = root.get("resolved_mass")
            if resolved is None and root.get("resolved_mass_num") is not None:
                resolved = f"{root['resolved_mass_num']}/{root['resolved_mass_den']}"
            if entry.get("excluded") and resolved is not None:
                known[p] = known.get(p, Fraction(0)) + as_fraction(resolved)
            elif status == "certified":
                if root.get("mass_num") is None or root.get("mass_den") is None:
                    raise AssertionError(f"certified root has no mass: {class_row['cell']} p={p}")
                known[p] = known.get(p, Fraction(0)) + Fraction(
                    int(root["mass_num"]), int(root["mass_den"])
                )
            elif resolved is not None:
                known[p] = known.get(p, Fraction(0)) + as_fraction(resolved)
            # A resolved root may still carry an unresolved residual bound;
            # accumulate it independently of the known mass.
            raw = root.get("residual_mass_upper_bound")
            if raw is None:
                raw = root.get("residual_mass_upper_bound_float")
            if raw is not None and entry.get("excluded"):
                residual[p] = residual.get(p, Fraction(0)) + as_fraction(raw)
    return known, residual, excluded


def rankdata(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for pos in order[i:j]:
            ranks[pos] = rank
        i = j
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("Spearman requires paired data")
    rx, ry = rankdata(x), rankdata(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    deny = math.sqrt(sum((b - my) ** 2 for b in ry))
    if denx == 0.0 or deny == 0.0:
        return float("nan")
    return num / (denx * deny)


def make_factors() -> tuple[list[dict], dict[tuple[int, tuple[int, int]], dict]]:
    source_rows = read_jsonl(BC_PATH)
    classes = [row for row in source_rows if row.get("type") == "class"]
    assert len(classes) == 293, len(classes)
    result: list[dict] = []
    exact_by_cell: dict[tuple[int, tuple[int, int]], dict] = {}

    for i, row in enumerate(classes, 1):
        key = cell_key(row["cell"])
        params = row["params"]
        ftable = {int(p): (int(num), int(den)) for p, num, den in row["f_p_table"]}
        assert len(ftable) == len(row["f_p_table"]), f"duplicate f_p_table p: {key}"
        known, residual, excluded = root_masses(row)

        # BAD_ROOT_LISTS is authoritative for this corrected artifact.  The
        # f_p_table entry for an excluded p can be blank, so it is audited
        # below but never used to manufacture clean mass.
        bad_ps = {int(entry["p"]) for entry in row["bad_root_lists"] if int(entry["p"]) <= CAP}
        assert set(p for p in ftable if p <= CAP) <= bad_ps, key
        for p, (num, den) in ftable.items():
            if p <= CAP:
                assert 0 <= num <= den and den > 0, (key, p, num, den)
                if p in excluded:
                    assert (num, den) == (0, 1), (key, p, (num, den))
                else:
                    assert Fraction(num, den) == known.get(p, Fraction(0)), (
                        key, p, (num, den), known.get(p, Fraction(0))
                    )

        all_window_ps = sorted(bad_ps | set(known) | set(residual))
        assert all(p <= CAP for p in all_window_ps)
        for p in all_window_ps:
            assert 0 <= known.get(p, Fraction(0)) <= 1, (key, p, known.get(p))
            assert 0 <= residual.get(p, Fraction(0)), (key, p, residual.get(p))
            assert known.get(p, Fraction(0)) + residual.get(p, Fraction(0)) <= 1, (
                key, p, known.get(p), residual.get(p)
            )
        upper_exact = math.prod(
            (1 - known.get(p, Fraction(0)) for p in all_window_ps),
            start=Fraction(1),
        )
        lower_exact = math.prod(
            (
                max(Fraction(0), 1 - known.get(p, Fraction(0)) - residual.get(p, Fraction(0)))
                for p in all_window_ps
            ),
            start=Fraction(1),
        )

        # The historical class-level field is intentionally not used here:
        # the corrected interval is built from each root's residual_p.
        top3 = sorted(
            ((p, known[p]) for p in known if p <= CAP and known[p] > 0),
            key=lambda item: (-item[1], item[0]),
        )[:3]
        out = {
            "cell": row["cell"],
            "family": row["family"],
            "a": params["a"],
            "eps": params["eps"],
            "f": params["f"],
            "q1": params["q1"],
            "N": params["N"],
            "partial": bool(row["partial"]),
            "factor_upper": float(upper_exact),
            "factor_lower": float(lower_exact),
            "factor_window": CAP,
            "top3_bad": [[p, float(mass)] for p, mass in top3],
        }
        assert set(out) == {
            "cell", "family", "a", "eps", "f", "q1", "N", "partial",
            "factor_upper", "factor_lower", "factor_window", "top3_bad",
        }
        result.append(out)
        exact_by_cell[key] = {
            "upper": upper_exact,
            "lower": lower_exact,
            "row": out,
            "top3": top3,
        }
        if i % 25 == 0:
            time.sleep(0.1)

    result.sort(key=lambda item: cell_key(item["cell"]))
    assert len(result) == 293
    with OUT_PATH.open("w") as handle:
        for item in result:
            handle.write(json.dumps(item, separators=(",", ":")) + "\n")

    # Verify the serialized artifact, not just the in-memory rows.
    reread = read_jsonl(OUT_PATH)
    assert len(reread) == 293
    for item in reread:
        exact = exact_by_cell[cell_key(item["cell"])]
        assert item["factor_window"] == CAP
        assert math.isclose(item["factor_upper"], float(exact["upper"]), rel_tol=0, abs_tol=2e-15)
        assert math.isclose(item["factor_lower"], float(exact["lower"]), rel_tol=0, abs_tol=2e-15)
        assert item["top3_bad"] == [[p, float(m)] for p, m in exact["top3"]]
    return result, exact_by_cell


def load_source(path: Path) -> list[dict]:
    return read_jsonl(path)


def closure_effort(records: list[dict], factors: dict[tuple[int, tuple[int, int]], dict]) -> tuple[list[dict], dict]:
    # Raw logs have different provenance fields.  Index each source once.
    caches: dict[str, dict[str, dict]] = {}
    for rec in records:
        source = rec.get("source")
        if source and source not in caches:
            rows = load_source(DATA / source)
            caches[source] = {
                "class": {cell_key(r["cell"]): r for r in rows if r.get("type") == "class" and r.get("cell")},
                "done": {cell_key(r["cell"]): r for r in rows if r.get("type") == "cell-done" and r.get("cell")},
            }

    l11 = [rec for rec in records if rec.get("family") == "L11"]
    assert len(l11) == 190, len(l11)
    joined: list[dict] = []
    source_counts: dict[str, int] = {}
    mode_counts: dict[str, int] = {}
    raw_fields: dict[str, dict] = {}
    for rec in l11:
        key = cell_key(rec["cell"])
        assert key in factors
        source = rec["source"]
        source_counts[source] = source_counts.get(source, 0) + 1
        cache = caches[source]
        class_row = cache["class"].get(key)
        done_row = cache["done"].get(key)
        is_alt = source.startswith("l13h_alt_l11")
        alt_tried = done_row.get("alt_classes_tried") if done_row else None
        if is_alt:
            raw_field = "cell-done.alt_classes_tried"
            raw_value = alt_tried
        elif class_row is not None and class_row.get("class_idx") is not None:
            raw_field = "class.class_idx"
            raw_value = class_row["class_idx"]
        else:
            raw_field = "no direct effort field"
            raw_value = None
        field_info = raw_fields.setdefault(source, {"field": raw_field, "values": []})
        assert field_info["field"] == raw_field, (source, field_info["field"], raw_field)
        if raw_value is not None:
            field_info["values"].append(int(raw_value))
        source_cls = class_row.get("class_idx") if class_row else (done_row.get("idx") if done_row else None)
        if is_alt:
            assert isinstance(alt_tried, int) and alt_tried >= 1, (source, rec["cell"], alt_tried)
            effort_cls = alt_tried - 1
            effort_mode = "alternate_classes_tried_minus_one"
        else:
            # Flat and deep-k use the frozen class itself; deep-k changes k,
            # not the class variant.  Deep-k has no direct cls field.
            effort_cls = 0
            effort_mode = "frozen_class_cls_zero"
        mode_counts[effort_mode] = mode_counts.get(effort_mode, 0) + 1
        joined.append({
            "cell": rec["cell"],
            "factor": factors[key]["row"]["factor_upper"],
            "effort_cls": effort_cls,
            "effort_mode": effort_mode,
            "source": source,
            "source_cls_index": source_cls,
            "alt_classes_tried": alt_tried,
            "k_zero": rec.get("k_zero"),
        })
        if len(joined) % 25 == 0:
            time.sleep(0.1)
    joined.sort(key=lambda row: cell_key(row["cell"]))
    assert len(joined) == 190
    return joined, {"source_counts": source_counts, "mode_counts": mode_counts, "raw_fields": raw_fields}


def render_report(rows: list[dict], joined: list[dict], source_info: dict) -> str:
    top10 = sorted(rows, key=lambda row: (row["factor_upper"], cell_key(row["cell"])))[:10]
    xs = [float(row["effort_cls"]) for row in joined]
    ys = [float(row["factor"]) for row in joined]
    rho = spearman(xs, ys)
    cls0 = [row["factor"] for row in joined if row["effort_cls"] == 0]
    cls1 = [row["factor"] for row in joined if row["effort_cls"] >= 1]
    assert cls0 and cls1

    lines: list[str] = []
    lines.append("L17 class-level risk quantification and closure-effort cross-validation")
    lines.append("factor_upper: product over authoritative BAD_ROOT_LISTS p<=10000 of (1-known_p)")
    lines.append("  known_p uses certified mass_num/mass_den and patched resolved_mass_num/resolved_mass_den;")
    lines.append("  f_p_table is audited against known non-excluded masses, but never trusted for excluded columns")
    lines.append(f"class rows: {len(rows)} (required 293; arithmetic reloaded and verified)")
    lines.append("factor_lower: product over p<=10000 of max(0, 1-known_p-residual_p), multiplicative per p")
    lines.append("  small_factor_padic and all p>10000 entries are excluded; factor_window is exactly 10000")
    lines.append("  historical excluded_mass_upper_bound is not used for endpoints; residual_p comes from each root")
    lines.append("  and neither residual_p nor any other field supplies a tail bound for primes p > 10000")
    lines.append("top3_bad: dominant known bad masses m_p, all with p <= 10000 (window-labeled)")
    lines.append("")
    lines.append("10 lowest factor classes (factor_upper; lower is shown when distinguishable)")
    lines.append("cell                 family  factor_upper       factor_lower       top3_bad (p,m_p)")
    lines.append("-------------------  ------  -----------------  -----------------  ---------------------------")
    for row in top10:
        cell = str(row["cell"])
        top = ", ".join(f"({p},{mass:.9g})" for p, mass in row["top3_bad"])
        lines.append(
            f"{cell:<19}  {row['family']:<6}  {row['factor_upper']:.15g}  "
            f"{row['factor_lower']:.15g}  {top}"
        )
    lines.append("")
    lines.append("Closure-effort cross-validation")
    lines.append("--------------------------------")
    lines.append("join: l13h_all_closures.json records joined to factor rows by cell")
    lines.append(f"L11 flat-family closure records used: {len(joined)} (required 190)")
    lines.append("raw effort fields observed:")
    for source, count in sorted(source_info["source_counts"].items()):
        lines.append(f"  {source}: {count}")
    lines.append("raw-field summary: flat class rows expose class_idx (79); alternate cell-done rows expose")
    lines.append("  alt_classes_tried (109); deep-k class rows expose neither direct effort field (2)")
    lines.append("raw per-source field values:")
    for source, info in sorted(source_info["raw_fields"].items()):
        values = info["values"]
        if values:
            lines.append(f"  {source}: {info['field']} range {min(values)}..{max(values)}")
        else:
            lines.append(f"  {source}: {info['field']}")
    lines.append("chosen effort axis: effort_cls, a 0-based closing class-variant index")
    lines.append("  flat/deep-k frozen-class records -> cls=0; alternate records -> alt_classes_tried-1")
    lines.append("  deep-k has no direct cls field in its raw rows; it is marked by mode but remains the frozen class")
    lines.append("effort-mode counts:")
    for mode, count in sorted(source_info["mode_counts"].items()):
        lines.append(f"  {mode}: {count}")
    lines.append("")
    lines.append(f"Spearman rho(effort_cls, factor_upper; p<=10000 window) = {rho:.9f}")
    lines.append("descriptive only: no p-value or inferential significance claim")
    lines.append(f"median factor_upper, cls=0: {statistics.median(cls0):.15g} (n={len(cls0)})")
    lines.append(f"median factor_upper, cls>=1: {statistics.median(cls1):.15g} (n={len(cls1)})")
    lines.append(f"effort counts: n_total={len(joined)}, cls=0={len(cls0)}, cls>=1={len(cls1)}")
    lines.append("")
    lines.append("VERDICT: PROVED arithmetic/window bookkeeping for all 293 factor rows; closure-effort association is descriptive EVIDENCE only.")
    return "\n".join(lines) + "\n"


def main() -> None:
    rows, exact = make_factors()
    with CLOSURES_PATH.open() as handle:
        closure_doc = json.load(handle)
    joined, source_info = closure_effort(closure_doc["records"], exact)
    report = render_report(rows, joined, source_info)
    REPORT.write_text(report)
    print(f"wrote {OUT_PATH} ({len(rows)} rows)")
    print(f"wrote {REPORT}")
    print(report.split("VERDICT:", 1)[-1].strip())


if __name__ == "__main__":
    main()
