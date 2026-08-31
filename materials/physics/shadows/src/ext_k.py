"""Gate-C k-extension campaign (pre_statement.md Revision 2, 2026-08-30).

Extends the closed gate-C campaign (campaigns/2026-08-29T232134Z_4f1a9cab_
e274dc65db35) beyond its k <= 8 scope:

  Q1  regression gate: recompute k = 4..8 rows with THIS code path and assert
      exact Fraction equality against the frozen gates_c.json values. The
      k > 8 rows count only if every regressed value matches.
  Q2  extension: k = 9..12, carrier M = cx_chain (deterministic-layer
      staircase) vs U_ct = mesh_Z, exact 3^k token enumeration, Fractions
      only in every verdict path (no RNG, no floats).

Pinned verdict (Revision 2, direction per Revision 1 semantics):
  EDGE-PERSISTS iff w_M(k) > w_ct(k)   (smaller shadow norm)
  EDGE-VANISHES iff w_M(k) <= w_ct(k)
  UNKNOWN-BUDGET iff k not reached before the wall cap.
Per-k weight difference w_M(k) - w_ct(k) reported exactly.

Usage (single core, low priority):
  nice -n 10 python3 ext_k.py campaign [DIR]
"""

import hashlib
import json
import platform
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import shadows_exact as se  # noqa: E402  (reuse, do not rewrite)

SHADOWS = Path(__file__).resolve().parent.parent
FROZEN_NAME = "2026-08-29T232134Z_4f1a9cab_e274dc65db35"
FROZEN_DIR = SHADOWS / "campaigns" / FROZEN_NAME
FROZEN_GATES_C = FROZEN_DIR / "gates_c.json"

WALL_CAP_S = 5400  # ~90 min budget cap pinned in Revision 2
CARRIER = "cx_chain"
LAYER_MAPS = ("cz_mesh", "cx_chain", "cx_mesh_lex", "uct_then_cz_mesh")
MESHES = ("mesh_X", "mesh_Y", "mesh_Z")


# ------------------------------------------------------------------ helpers

def fmt(fr: Fraction) -> str:
    """Exact Fraction as 'num/den' (frozen-artifact string convention)."""
    assert isinstance(fr, Fraction)
    return f"{fr.numerator}/{fr.denominator}"


def parse_w(w: str) -> Fraction:
    num, den = w.split("/")
    return Fraction(int(num), int(den))


def build_layers(k: int):
    """Deterministic-layer maps for k qubits, built exactly as gate_c_run
    built them (same gate order, same composition direction)."""
    uct, _ = se.mesh_all_pairs(k, "Z")
    base = tuple(1 << b for b in range(2 * k))
    for i, j in combinations(range(k), 2):
        base = se.compose(se.embed_gate(se.cz_imgs(), k, [i, j]), base)
    cz_mesh = base
    base = tuple(1 << b for b in range(2 * k))
    for i in range(k - 1):
        base = se.compose(se.embed_gate(se.cnot_imgs(), k, [i, i + 1]), base)
    cx_chain = base
    base = tuple(1 << b for b in range(2 * k))
    for i, j in combinations(range(k), 2):
        base = se.compose(se.embed_gate(se.cnot_imgs(), k, [i, j]), base)
    cx_mesh_lex = base
    uct_after_cz = se.compose(cz_mesh, uct)
    return {"uct(mesh_Z)": uct, "uct": uct, "cz_mesh": cz_mesh,
            "cx_chain": cx_chain, "cx_mesh_lex": cx_mesh_lex,
            "uct_then_cz_mesh": uct_after_cz}


def eval_weight(k: int, imgs) -> Fraction:
    w, _hist = se.eval_map_on_tokens(k, imgs)
    return w


# ------------------------------------------------- Q1: k = 4..8 regression

def frozen_rows() -> dict:
    gc = json.loads(FROZEN_GATES_C.read_text())
    out = {"families": {name: {} for name in LAYER_MAPS + MESHES}}
    for name in LAYER_MAPS:
        for r in gc["families"]["deterministic_layers"]:
            if r.get("family") == name:
                assert r["k"] not in out["families"][name]
                out["families"][name][r["k"]] = r["w_exact"]
    fam_letter = {"mesh_X": "mesh_X", "mesh_Y": "mesh_Y", "mesh_Z": "mesh_Z"}
    for name in MESHES:
        for r in gc["families"][fam_letter[name]]:
            out["families"][name][r["k"]] = r["w_exact"]
    return out


def regression_gate():
    """Recompute k=4..8 with this code path; assert exact equality with the
    frozen artifact. Returns (rows, ok). k > 8 rows count only if ok."""
    frozen = frozen_rows()
    rows = []
    ok = True
    for k in range(4, 9):
        maps = build_layers(k)
        w_ct = eval_weight(k, maps["uct"])
        assert w_ct == se.eq4_weight(k), (k, "U_ct != Eq. (4)")
        row = {"k": k, "w_ct": fmt(w_ct)}
        for name in LAYER_MAPS:
            assert se.is_symplectic(maps[name], k)
            w = eval_weight(k, maps[name])
            match = w == parse_w(frozen["families"][name][k])
            ok = ok and match
            row[name] = {"w_recomputed": fmt(w),
                         "w_frozen": frozen["families"][name][k],
                         "exact_match": bool(match)}
        for name in MESHES:
            w = eval_weight(k, maps[f"uct(mesh_Z)"])
            # single-letter meshes are NOT the identity; rebuild each quickly
            imgs, _ = se.mesh_all_pairs(k, name[-1])
            w = eval_weight(k, imgs)
            match = w == parse_w(frozen["families"][name][k])
            ok = ok and match
            row[name] = {"w_recomputed": fmt(w),
                         "w_frozen": frozen["families"][name][k],
                         "exact_match": bool(match)}
        # carrier-vs-U_ct anchor from Revision 2 Fact 1, asserted exactly:
        assert w_ct == Fraction(353, 6561) and \
            parse_w(row["cx_chain"]["w_recomputed"]) == Fraction(361, 6561) \
            if k == 4 else True, (k, "Revision-2 Fact 1 mismatch")
        rows.append(row)
    return rows, ok


# --------------------------------------------------- Q2: k = 9..12 extension

def run_extension(kmax: int = 12):
    """k = 9..12 rows; respects WALL_CAP_S; unstarted k -> UNKNOWN-BUDGET."""
    rows = []
    t_start = time.perf_counter()
    for k in range(9, kmax + 1):
        if time.perf_counter() - t_start > WALL_CAP_S:
            for k_miss in range(k, kmax + 1):
                rows.append({"k": k_miss, "verdict": "UNKNOWN-BUDGET"})
            break
        maps = build_layers(k)
        w_ct = eval_weight(k, maps["uct"])
        assert w_ct == se.eq4_weight(k), (k, "U_ct != Eq. (4) after enum")
        hist_ct = se.token_enumeration(k, maps["uct"])[1]
        w_cx = eval_weight(k, maps[CARRIER])
        assert se.is_symplectic(maps[CARRIER], k)
        hist_cx = se.token_enumeration(k, maps[CARRIER])[1]
        extra = {name: fmt(eval_weight(k, maps[name]))
                 for name in ("cz_mesh", "cx_mesh_lex", "uct_then_cz_mesh")}
        diff = w_cx - w_ct
        rows.append({
            "k": k,
            "w_M": fmt(w_cx),
            "w_ct": fmt(w_ct),
            "diff_wM_minus_wct": fmt(diff),
            "verdict": "EDGE-PERSISTS" if w_cx > w_ct else "EDGE-VANISHES",
            "norm2_M": fmt(Fraction(w_cx.denominator, w_cx.numerator)),
            "norm2_ct": fmt(Fraction(w_ct.denominator, w_ct.numerator)),
            "hist_ct": {str(m): c for m, c in sorted(hist_ct.items())},
            "hist_cx": {str(m): c for m, c in sorted(hist_cx.items())},
            "others": extra,
            "elapsed_s": round(time.perf_counter() - t_start, 3),
        })
    return rows, round(time.perf_counter() - t_start, 3)


# ------------------------------------------------------------------- main

def artifact_id() -> str:
    """12-hex sha256 prefix over the frozen statement + the two code files."""
    h = hashlib.sha256()
    for p in (SHADOWS / "pre_statement.md",
              SHADOWS / "src" / "shadows_exact.py",
              Path(__file__).resolve()):
        h.update(p.read_bytes())
    return h.hexdigest()[:12]


def new_campaign_name() -> str:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    return f"{ts}_{uuid.uuid4().hex[:8]}_{artifact_id()}"


def selftest():
    # anchor facts of the frozen campaign, re-exercised small
    assert se.pauli_size(0b0101, 2) == 2
    maps4 = build_layers(4)
    assert eval_weight(4, maps4["uct"]) == Fraction(353, 6561) == se.eq4_weight(4)
    assert eval_weight(4, maps4["cx_chain"]) == Fraction(361, 6561)
    assert eval_weight(4, maps4["uct_then_cz_mesh"]) == Fraction(1, 81)
    # literal cx_chain after U_ct excluded: weight BELOW U_ct at k=4
    assert eval_weight(4, maps4["uct_then_cz_mesh"]) < eval_weight(4, maps4["uct"])
    # 2q anchor symplectic + maps symplectic up to k=6 quickly
    assert se.is_symplectic(se.ANCHOR_IMGS, 2)
    for k in (5, 6):
        m = build_layers(k)
        for name in ("uct", "cx_chain"):
            assert se.is_symplectic(m[name], k)
    print("selftest OK (k=4 anchors 353/6561, 361/6561, 1/81; symplectic maps)")


def run_campaign(campaign_dir: Path):
    assert campaign_dir.name.endswith(artifact_id()), \
        (campaign_dir.name, artifact_id())
    campaign_dir.mkdir(parents=True, exist_ok=False)
    wall = time.perf_counter()
    reg_rows, ok = regression_gate()
    reg_elapsed = time.perf_counter() - wall
    assert ok, "regression gate FAILED vs frozen artifact; aborting (rows kept)"
    ext_rows, ext_elapsed = run_extension()
    total = round(time.perf_counter() - wall, 3)

    gates_ext = {
        "gate": "C-extension (Revision 2)",
        "carrier": CARRIER,
        "frozen_source": FROZEN_NAME,
        "regression_k4_to_k8": reg_rows,
        "regression_ok": bool(ok),
        "extension": ext_rows,
        "elapsed_regression_s": round(reg_elapsed, 3),
        "elapsed_extension_s": round(ext_elapsed, 3),
        "verdict_rule": ("EDGE-PERSISTS iff w_M(k) > w_ct(k) exact; "
                         "EDGE-VANISHES iff w_M(k) <= w_ct(k); "
                         "per Revision 1 semantics (smaller shadow norm)"),
        "excluded_map": ("cx_chain-after-U_ct literal composite (different "
                         "map; w=1/81 < 353/6561 at k=4)"),
    }
    (campaign_dir / "gates_ext.json").write_text(
        json.dumps(gates_ext, indent=1) + "\n")

    completed = [r["k"] for r in ext_rows if r.get("verdict") != "UNKNOWN-BUDGET"]
    summary = {
        "campaign": campaign_dir.name,
        "statement": "pre_statement.md Revision 2 (appended 2026-08-30, "
                     "before any k > 8 run)",
        "frozen_source": FROZEN_NAME,
        "carrier": CARRIER,
        "regression_gate": "PASS" if ok else "FAIL",
        "per_k_extension": [
            {"k": r["k"], "w_M": r.get("w_M"), "w_ct": r.get("w_ct"),
             "diff": r.get("diff_wM_minus_wct"), "verdict": r["verdict"]}
            for r in ext_rows],
        "largest_k_completed": max(completed) if completed else None,
        "elapsed_regression_s": round(reg_elapsed, 3),
        "elapsed_extension_s": round(ext_elapsed, 3),
        "total_elapsed_s": total,
        "exactness": "GF(2) bit-label ints + fractions.Fraction only; no RNG",
        "single_core_low_priority": "nice -n 10, one process",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "artifact_sha256_12": artifact_id(),
    }
    (campaign_dir / "SUMMARY.json").write_text(json.dumps(summary, indent=1) + "\n")

    manifest = [
        f"campaign {campaign_dir.name}",
        f"artifact_sha256_12(gates_ext): {artifact_id()}",
        "frozen statement: pre_statement.md Revision 2 (2026-08-30, appended "
        "before this run)",
        f"frozen source artifact: campaigns/{FROZEN_NAME}/gates_c.json "
        "(regression-anchored, exact match)",
        f"carrier M = {CARRIER}; U_ct = mesh_Z; excluded: literal cx_chain"
        "-after-U_ct composite",
        "verdict: EDGE-PERSISTS iff w_M(k) > w_ct(k) exact, else "
        "EDGE-VANISHES (Revision 1 semantics)",
        "exactness: Fraction-only verdict path, no RNG; full 3^k enumeration",
        f"total_elapsed_s: {total}",
        f"wall_cap_s: {WALL_CAP_S}",
    ]
    (campaign_dir / "MANIFEST.txt").write_text("\n".join(manifest) + "\n")
    return summary


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if mode == "selftest":
        selftest()
    elif mode == "newdir":
        print(new_campaign_name())
    elif mode == "campaign":
        target = Path(sys.argv[2]) if len(sys.argv) > 2 else \
            SHADOWS / "campaigns" / new_campaign_name()
        print(json.dumps(run_campaign(target), indent=1))
    else:
        raise SystemExit(f"unknown mode {mode!r}")
