"""EXP-052: ideal-power (quotient-action) classification of the demote law.

For each catalogue BB parent P with X-check block H_X = [A B] (lattice
(ell, m), ring R = GF(2)[x,y]/(x^ell-1, y^m-1)):

  M = ker H_X / S_Z            (k_P-dimensional quotient; same coordinates as
                                EXP-047's ``quotient_and_projector``)
  I = L_pre = nullspace(H_X.T) (the left-null ideal, dim k_P/2 by Lemma 2)

The demote-fixed set S = {y in M : y in I y} is the stable kernel of the
ideal-power chain

    F_0 = M,   F_{r+1} = I F_r

because R is finite commutative, hence Artinian: write R = prod_beta R_beta
with LOCAL Artinian factors R_beta (note: the coarse (p^a, q^b) primary
blocks are NOT in general local -- e.g. GF(4) ⊗ GF(4) = GF(4) x GF(4) on the
6x6 lattice -- so idempotent analysis is done only through the F-chain, which
needs no factor decomposition).  Then F_infty = e_I M where e_I is the
idempotent of I's full-factor support, and

    y in I y  <=>  y is supported on I-full factors  <=>  y = e_I y.

(PROOF INDEPENDENTLY AUDITED 2026-08-21 by codex CLI read-only review: all
algebraic steps confirmed; three presentation defects fixed in
notes/theorem_je2_demote_trichotomy.md: local-unit phrasing in the local
step, the nonzero-class fraction convention, and the descent/stop-rule
wording.)

Classification per parent:
    dim S = k_P  -> immune        (demote fraction 0)
    dim S = 0    -> demote-full   (demote fraction 1)
    else         -> mixed         (demote fraction 1 - (2^{dim S}-1)/(2^{k_P}-1))

Cross-checks (run mode persists per-parent JSON; assemble verifies):
  * immune set == EXP-047 immunity_exact list (10 parents),
  * fraction_pred == EXP-050 fraction on all 196 enumerated parents,
  * the six k_P>20 parents (enumeration-infeasible) are classified here for
    the first time.

This replaces enumeration over 2^{k_P}-1 quotient classes by O(chain depth)
GF(2) linear algebra and extends the machine census from 196/202 to all 202
parents.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC27 = importlib.util.spec_from_file_location("exp027", ROOT / "experiments" / "exp027_delta_audit.py")
assert _SPEC27 and _SPEC27.loader
E27 = importlib.util.module_from_spec(_SPEC27)
sys.modules[_SPEC27.name] = E27
_SPEC27.loader.exec_module(E27)

_SPEC39 = importlib.util.spec_from_file_location("exp039", ROOT / "experiments" / "exp039_nogo_module.py")
assert _SPEC39 and _SPEC39.loader
E39 = importlib.util.module_from_spec(_SPEC39)
sys.modules[_SPEC39.name] = E39
_SPEC39.loader.exec_module(E39)

_SPEC47 = importlib.util.spec_from_file_location("exp047", ROOT / "experiments" / "exp047_exact_demotion_decision.py")
assert _SPEC47 and _SPEC47.loader
E47 = importlib.util.module_from_spec(_SPEC47)
sys.modules[_SPEC47.name] = E47
_SPEC47.loader.exec_module(E47)

from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np  # noqa: E402

SCHEMA = "exp052-ideal-power-trichotomy-v1"
STATE_DIR = ROOT / "results" / "partial_runs" / "exp052"
OUT = ROOT / "results" / "processed" / "exp052_ideal_power_trichotomy.json"
MAX_CHAIN_STEPS = 200


def classify_case(s: int, k_P: int) -> tuple[str, float]:
    """Map the stable dimension dim S = F_infty to the trichotomy case and the
    predicted demote fraction 1 - (2^{s}-1)/(2^{k_P}-1)."""
    if s == k_P:
        return "immune", 0.0
    if s == 0:
        return "demote_full", 1.0
    return "mixed", 1.0 - ((1 << s) - 1) / ((1 << k_P) - 1)


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def atomic_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=1, sort_keys=True))
    tmp.replace(path)


def parent_chain(fp: str, entry: dict, row: dict) -> dict:
    """Compute the F-chain for one parent and classify it."""
    t0 = time.time()
    ell, m = int(entry["ell"]), int(entry["m"])
    dim = ell * m
    _, HX, HZ = E27.parent_matrices(row)
    k_P = int(entry["k_parent"])

    Kb = nullspace_np(HX)
    Q, R, kq = E47.quotient_and_projector(HZ, Kb)
    if kq != k_P:
        raise RuntimeError(f"{entry['members'][0]['label']}: quotient dim {kq} != k_P {k_P}")
    Rk = R[:, :k_P]
    Lpre = nullspace_np(HX.T)
    if 2 * Lpre.shape[0] != k_P:
        raise RuntimeError(f"{entry['members'][0]['label']}: dim L_pre {Lpre.shape[0]} != k_P/2")

    # F_0 = M in coordinates: identity basis
    basis = np.eye(k_P, dtype=np.uint8)
    chain = [k_P]
    while True:
        rows = []
        for b in basis:
            z = (b @ Q) % 2
            sz = E47.shift_matrix_of(z, ell, m)
            rows.append((Lpre @ sz % 2) @ Rk % 2)
        nxt, _ = rref_np(np.vstack(rows)) if rows else (np.zeros((0, k_P), np.uint8), [])
        r_next = rank_np(nxt)
        r_cur = rank_np(basis)
        # Containment F_{r+1} <= F_r is mathematically automatic (ideal powers
        # nest); the joint-rank assertion below is a regression tripwire: a
        # wrongly coded quotient action that violates the monotone law dies
        # loudly here instead of stabilizing at an equal-rank wrong subspace.
        if nxt.size:
            joint = rank_np(np.vstack([basis, nxt]))
        else:
            joint = r_cur
        if r_next > r_cur or joint != r_cur:
            raise RuntimeError(
                f"{entry['members'][0]['label']}: F-chain violated monotonicity "
                f"(rank {r_cur} -> {r_next}, joint {joint})"
            )
        chain.append(int(r_next))
        if r_next == r_cur:
            break  # automatic containment + equal rank => F_{r+1} = F_r
        basis = nxt
        if len(chain) > MAX_CHAIN_STEPS:
            raise RuntimeError(f"{entry['members'][0]['label']}: F-chain did not stabilize")

    s = chain[-1]
    case, fraction_pred = classify_case(s, k_P)

    return {
        "schema": SCHEMA,
        "fingerprint": fp,
        "label": entry["members"][0]["label"],
        "ell": ell,
        "m": m,
        "n": 2 * dim,
        "k_parent": k_P,
        "dim_L_pre": int(Lpre.shape[0]),
        "chain": chain,
        "chain_steps": len(chain) - 1,
        "dim_S": int(s),
        "case": case,
        "fraction_pred": fraction_pred,
        "wall_s": round(time.time() - t0, 3),
        "utc": utc_now(),
    }


def run(args: argparse.Namespace) -> int:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    only = set(args.only or [])
    done = 0
    for fp, entry in sorted(parents.items(), key=lambda kv: (kv[1]["n"], kv[1]["k_parent"], kv[0])):
        if only and fp not in only and entry["members"][0]["label"] not in only:
            continue
        out = STATE_DIR / f"chain_{fp}.json"
        if out.exists() and not args.force:
            done += 1
            continue
        rec = parent_chain(fp, entry, rows[entry["members"][0]["catalogue_index"]])
        atomic_write_json(out, rec)
        done += 1
        print(f"[{len(parents)}] {rec['label']}: k_P={rec['k_parent']} chain={rec['chain']} case={rec['case']}",
              flush=True)
    print(f"parents processed or cached: {done}")
    return 0


def assemble(args: argparse.Namespace) -> int:
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    recs = {}
    for fp, entry in parents.items():
        p = STATE_DIR / f"chain_{fp}.json"
        if not p.exists():
            raise RuntimeError(f"missing chain record for {fp} ({entry['members'][0]['label']})")
        recs[fp] = json.loads(p.read_text())

    cases = {}
    for rec in recs.values():
        cases[rec["case"]] = cases.get(rec["case"], 0) + 1

    d47 = json.loads((ROOT / "results/processed/exp047_exact_demotion_decision.json").read_text())
    immune47 = {fp for fp, _ in d47["immunity_list"]}
    pred_immune = {fp for fp, r in recs.items() if r["case"] == "immune"}
    # long vs 12-hex fingerprints
    long_of_12 = {fp[:12]: fp for fp in parents}
    immune47_long = {long_of_12[x] for x in immune47}
    immune_mismatch = sorted((pred_immune ^ immune47_long))

    fr_files = sorted(STATE_DIR.parent.glob("exp050/fraction_*.json"))
    if len(fr_files) != 196:
        raise RuntimeError(f"expected 196 exp050 records, found {len(fr_files)}")
    fraction_mismatches = []
    crosschecked = 0
    for f in fr_files:
        d = json.loads(f.read_text())
        fp = d["fingerprint"]
        if fp not in recs:
            continue
        crosschecked += 1
        # Exact integer comparison: classes_demote == 2^k_P - 2^{dim S}
        k = recs[fp]["k_parent"]
        s = recs[fp]["dim_S"]
        expected = (1 << k) - (1 << s)
        if int(d["classes_demote"]) != expected or int(d["classes_total"]) != (1 << k) - 1:
            fraction_mismatches.append({
                "fingerprint": fp, "label": d["label"],
                "exp050": [int(d["classes_demote"]), int(d["classes_total"])],
                "pred": [expected, (1 << k) - 1],
            })
    if crosschecked != 196:
        raise RuntimeError(f"expected 196 crosschecked parents, got {crosschecked}")

    heavy = sorted([r["label"] for r in recs.values() if r["k_parent"] > 20])
    heavy_cases = {r["label"]: r["case"] for r in recs.values() if r["k_parent"] > 20}
    max_chain = max(r["chain_steps"] for r in recs.values())

    out = {
        "schema": SCHEMA + "-assembled",
        "experiment": "exp052_ideal_power_trichotomy",
        "utc": utc_now(),
        "parents_total": len(recs),
        "cases": cases,
        "immune_labels": sorted(r["label"] for r in recs.values() if r["case"] == "immune"),
        "mixed_labels": sorted(r["label"] for r in recs.values() if r["case"] == "mixed"),
        "heavy_gt20": heavy_cases,
        "max_chain_steps": max_chain,
        "immune_mismatch_vs_exp047": immune_mismatch,
        "crosschecked_fractions": crosschecked,
        "fraction_mismatches_vs_exp050": fraction_mismatches,
        "records": recs,
    }
    atomic_write_json(OUT, out)
    print(json.dumps({k: out[k] for k in (
        "parents_total", "cases", "mixed_labels", "heavy_gt20", "max_chain_steps",
        "immune_mismatch_vs_exp047", "crosschecked_fractions")}, indent=1))
    print(f"fraction_mismatches_vs_exp050: {len(fraction_mismatches)}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("run")
    pr.add_argument("--only", nargs="*", default=None)
    pr.add_argument("--force", action="store_true")
    pr.set_defaults(fn=run)
    pa = sub.add_parser("assemble")
    pa.set_defaults(fn=assemble)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
