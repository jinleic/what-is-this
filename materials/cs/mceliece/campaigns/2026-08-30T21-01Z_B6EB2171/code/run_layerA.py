"""Gate B Layer A runner — sampled seed stream over the registered cells.

Cells (ADDENDUM 1; base pre_statement): m in {6..11}, small t per m:
  m=6: t in {2,3}; m=7: t in {2,3}; m=8: t in {2,3}; m=9: t in {2,3};
  m=10: t = 3; m=11: t = 3.   (t=1 at full support is build-impossible:
  the degree-1 irreducible G has its unique root IN the full support.)

Seeds: 25 per cell, stream 1000..1024 per cell (registered). A build whose
guards fail is logged with its failure and the stream ADVANCES to the next
seed (Gate-A pre_statement §6 rule); the abort is recorded in the results.

Checkpointing after EVERY instance (two predecessors lost a full run with
nothing on disk): state.json rewritten after each build; append-only log
per cell. Verdict per instance: engine_verdict (exact; pointscan route with
witness or pointscan_complete / all_squares). Every DEGENERATE verdict is
re-verified by the exact brute Wronskian over ALL pairs before being
recorded as an event.

nice-10, single process, threads pinned to 1 (repo resource policy).
"""
from __future__ import annotations
import json, os, sys, time, hashlib

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from gfield import GF
from fastfield import EField
from instance import Instance
from engine import engine_verdict, cascade_verdict

OUT = os.path.dirname(HERE)          # scratch/gateB
STATE = os.path.join(OUT, "layerA_state.json")
LOG = os.path.join(OUT, "layerA_log.jsonl")

CELLS = []
for m, ts in [(6, (2, 3)), (7, (2, 3)), (8, (2, 3)), (9, (2, 3)), (10, (3,)), (11, (3,))]:
    for t in ts:
        CELLS.append({"m": m, "t": t, "n": 1 << m})
SEEDS_PER_CELL = 25
STREAM_BASE = 1000


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_state():
    if os.path.exists(STATE):
        return json.load(open(STATE))
    return {"cells": {}, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


def save_state(st):
    tmp = STATE + ".tmp"
    json.dump(st, open(tmp, "w"), indent=1)
    os.replace(tmp, STATE)


def append_log(rec):
    with open(LOG, "a") as f:
        f.write(json.dumps(rec, default=str) + "\n")


def run_cell(cell, st):
    key = f"m{cell['m']}_t{cell['t']}"
    if key in st["cells"] and st["cells"][key].get("done"):
        print(f"[skip] {key} done", flush=True)
        return
    m, n, t = cell["m"], cell["n"], cell["t"]
    gf = GF(m)
    ef = EField(m)
    cstate = st["cells"].setdefault(key, {"records": [], "next": 0, "done": False,
                                          "aborts": []})
    t_cell = time.time()
    while cstate["next"] < SEEDS_PER_CELL and not cstate["done"]:
        sidx = cstate["next"]
        seed = STREAM_BASE + sidx
        t0 = time.time()
        rec = {"cell": key, "m": m, "n": n, "t": t, "seed": seed,
               "seed_index": sidx}
        try:
            inst = Instance(m, n, t, seed)
            guards = inst.guards
            rec["guards_ok"] = all(guards.get(x, False) for x in
                                   ("alpha_identity", "gamma_k", "eps_gcd_const",
                                    "eps_maxdeg_is_D"))
            if not rec["guards_ok"]:
                rec["verdict"] = "GUARD_FAIL"
                rec["guard_detail"] = {kk: vv for kk, vv in guards.items()
                                       if kk != "degs"}
                cstate["aborts"].append({"seed": seed,
                                         "detail": rec["guard_detail"]})
                # advance stream (registered rule)
                cstate["next"] += 1
                cstate["records"].append(rec)
                append_log(rec)
                save_state(st)
                print(f"[{key}] seed {seed}: GUARD_FAIL (stream advances)",
                      flush=True)
                continue
            ev = engine_verdict(inst.F, inst.Y, np.asarray(inst.lam, np.uint16),
                                inst.support, ef, gf)
            rec["engine"] = {kk: ev[kk] for kk in ev if kk != "beta_at_points"}
            if ev["verdict"] == "DEGENERATE":
                # EVENT: re-verify with exact brute Wronskian over ALL pairs
                cv = cascade_verdict(gf, [gf.ptrim(list(f)) for f in inst.F])
                rec["brute_recheck"] = cv
                if cv["verdict"] != "DEGENERATE":
                    rec["verdict"] = "INSTRUMENT_DISAGREEMENT"
                    raise RuntimeError(f"instrument disagreement at {key} seed {seed}: {ev} vs {cv}")
                rec["verdict"] = "DEGENERATE"
                rec["G"] = inst.G
                rec["g_tries"] = inst.g_tries
                rec["support"] = inst.support
                rec["lam"] = [int(x) for x in inst.lam]
                rec["F_degrees"] = [gf.pdeg(f) for f in inst.F]
            elif ev["verdict"] == "NONDEGENERATE":
                rec["verdict"] = "NONDEGENERATE"
            else:
                rec["verdict"] = ev["verdict"]
            cstate["next"] += 1
            cstate["records"].append(rec)
            append_log(rec)
            save_state(st)
            print(f"[{key}] seed {seed}: {rec['verdict']} "
                  f"({ev.get('route', '-')}) [{time.time()-t0:.1f}s]",
                  flush=True)
        except Exception as e:
            rec["verdict"] = "BUILD_ERROR"
            rec["error"] = repr(e)[:300]
            cstate["records"].append(rec)
            append_log(rec)
            save_state(st)
            print(f"[{key}] seed {seed}: BUILD_ERROR {e!r}", flush=True)
            cstate["next"] += 1
            continue
        # registered stop: on ANY degenerate event, halt and escalate
        if rec["verdict"] == "DEGENERATE":
            print(f"[{key}] DEGENERATE EVENT at seed {seed} — HALT for escalation",
                  flush=True)
            cstate["done"] = True
            st["halt_event"] = rec
            save_state(st)
            return
    cstate["done"] = cstate["next"] >= SEEDS_PER_CELL
    cstate["elapsed_s"] = round(time.time() - t_cell, 1)
    save_state(st)


def main():
    st = load_state()
    for cell in CELLS:
        if st.get("halt_event"):
            break
        run_cell(cell, st)
    print(json.dumps({"status": "halt" if st.get("halt_event") else "complete",
                      "cells": {k: {"n_records": len(v.get("records", [])),
                                    "done": v.get("done")}
                                for k, v in st["cells"].items()}},
                     indent=1))


if __name__ == "__main__":
    main()
