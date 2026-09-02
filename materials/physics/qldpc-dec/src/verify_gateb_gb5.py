"""Independent acceptance check for the GB5a, GB6, and GB7 campaigns.

Re-derives every reported number from the frozen raw artifacts rather than
trusting the summaries: shot-file hashes, prediction hashes, failure masks,
shard stitching, independent-arm ratio CIs, the three-way band decisions,
paired discordance, exact McNemar decisions, prefix continuity, and GB7's
paired factor intervals. Prints one PASS line per verified claim.

usage: verify_gateb_gb5.py <ladder_dir> [<mechanism_dir> [<gb7_dir>]]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

TARGET = Path(__file__).resolve().parents[1]
SRC = TARGET / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gateb_gb5a_ladder import (  # noqa: E402
    RUNGS, failure_mask_from_predictions, ladder_verdict, mcnemar_exact,
    ratio_verdict, read_failmask, sha256_file,
)
from gateb_gb7_paired import (  # noqa: E402
    BEAM_BINARY, PREFIX_CAMPAIGN, PREFIX_SHOTS, compare_prediction_prefix,
    compare_shot_prefix, paired_cells, rung_decision,
)
from gateb_gb8_retarget import (  # noqa: E402
    ARTIFACT as GB8_ARTIFACT, PAPER as GB8_PAPER, retarget_all,
)
from qldpc_dec.seeds import derive_seed  # noqa: E402

BASE_SEED = 20260829
failures: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    if ok:
        print(f"PASS  {label}" + (f" | {detail}" if detail else ""))
    else:
        failures.append(label)
        print(f"FAIL  {label}" + (f" | {detail}" if detail else ""))


def check_frozen_hashes(prefix: str, cdir: Path) -> None:
    """Post-close: every file in sha256s.txt must still hash identically.

    Pre-close there is nothing to compare against; that is reported as an
    explicit SKIP (never a PASS) so a pre-close run cannot masquerade as
    the final post-close acceptance.
    """
    listing = cdir / "sha256s.txt"
    if not listing.is_file():
        print(f"SKIP  {prefix}.frozen_hashes | not frozen yet; re-run after campaign.py close")
        return
    mismatched = []
    count = 0
    for line in listing.read_text().splitlines():
        digest, _sep, relative = line.partition("  ")
        count += 1
        if sha256_file(cdir / relative) != digest:
            mismatched.append(relative)
    check(f"{prefix}.frozen_hashes", not mismatched,
          f"{count} files re-hashed" + (f"; MISMATCH {mismatched}" if mismatched else ""))
    check(f"{prefix}.terminal_verdict", (cdir / "status.json").is_file(),
          json.loads((cdir / "status.json").read_text()).get("verdict", "?")
          if (cdir / "status.json").is_file() else "status.json missing")


def verify_ladder(cdir: Path) -> dict:
    manifest = json.loads((cdir / "manifest.json").read_text())
    summary = json.loads((cdir / "summary.json").read_text())
    inventory = json.loads((cdir / "sample_inventory.json").read_text())
    shots = int(summary["shots_per_arm"])

    check("ladder.campaign_closed", (cdir / "inventory.json").exists(), cdir.name)
    check_frozen_hashes("ladder", cdir)
    check("ladder.prior_channel_exact",
          manifest["bposd_prior_check"]["configured_channel_exact"] is True,
          f"{manifest['bposd_prior_check']['columns']} columns, "
          f"{manifest['bposd_prior_check']['distinct_probabilities']} distinct")

    eq = manifest["beam_equivalence_evidence"]
    check("ladder.equivalence_gate", len(eq) == 6 and all(e["mismatched_shots"] == 0 for e in eq),
          f"{len(eq)} reports, all zero-mismatch")
    for e in eq:
        p = TARGET / e["path"]
        check(f"ladder.equivalence_report_hash[{e['rung']}@p{e['p']:g}]",
              sha256_file(p) == e["sha256"], e["sha256"][:12])

    for arm, meta in inventory["arms"].items():
        path = cdir / str(meta["path"])
        check(f"ladder.sample_hash[{arm}]", sha256_file(path) == meta["sha256"],
              f"{meta['shots']} shots, {meta['sha256'][:12]}")

    bposd = summary["bposd"]
    shard_files = sorted(cdir.glob("bposd_shard_*.json"))
    shard_data = [json.loads(p.read_text()) for p in shard_files]
    total = sum(int(s["shots"]) for s in shard_data)
    fsum = sum(int(s["failures"]) for s in shard_data)
    check("ladder.bposd_shard_coverage", total == shots, f"{len(shard_data)} shards, {total} shots")
    check("ladder.bposd_failure_sum", fsum == int(bposd["failures"]), f"{fsum} failures")
    stitched = np.concatenate([read_failmask(cdir / s["failure_mask"], int(s["shots"]))
                               for s in sorted(shard_data, key=lambda s: int(s["shard"]))])
    check("ladder.bposd_mask_consistency",
          stitched.size == shots and int(stitched.sum()) == fsum,
          f"{int(stitched.sum())} set bits")

    masks = {"bp30+osd": stitched}
    for rung in RUNGS:
        arm = summary["rung_arms"][rung]
        pred = cdir / f"{rung}_predictions.bin"
        check(f"ladder.predictions_hash[{rung}]",
              sha256_file(pred) == arm["decode"]["predictions_sha256"], arm["decode"]["predictions_sha256"][:12])
        cfg = arm["decode"]["config_line"]
        want = (f"config: beam_width={RUNGS[rung]['beam_width']} "
                f"initial_iters={RUNGS[rung]['initial_iters']} "
                f"iters_per_round={RUNGS[rung]['iters_per_round']} "
                f"max_rounds={RUNGS[rung]['max_rounds']} num_results=1")
        check(f"ladder.decoder_config[{rung}]", cfg == want, cfg[8:])
        recomputed, _po = failure_mask_from_predictions(cdir / f"shots_{rung}.bin", pred)
        check(f"ladder.failures_recomputed[{rung}]",
              int(recomputed.sum()) == int(arm["failures"]),
              f"{int(recomputed.sum())} of {shots}")

    rung_verdicts = {}
    for rung in RUNGS:
        arm = summary["rung_arms"][rung]
        stored = summary["rung_verdicts"][rung]
        seed = derive_seed(BASE_SEED, "gateB-gb5a-ratio-bootstrap", rung, shots)
        redone = ratio_verdict(arm, bposd, RUNGS[rung]["band"], seed)
        rung_verdicts[rung] = redone
        same = (redone["outcome"] == stored["outcome"]
                and redone.get("beam_over_bposd_ratio") == stored.get("beam_over_bposd_ratio")
                and redone["ratio_ci95"] == stored["ratio_ci95"])
        check(f"ladder.verdict_reproducible[{rung}]", same,
              f"{redone['outcome']} ratio={redone.get('beam_over_bposd_ratio')} CI={redone['ratio_ci95']}")

    redone_ladder = ladder_verdict(rung_verdicts)
    check("ladder.ladder_verdict_reproducible",
          redone_ladder["ladder_outcome"] == summary["ladder_verdict"]["ladder_outcome"],
          redone_ladder["ladder_outcome"])

    paired = summary["paired_instrument"]
    for rung in RUNGS:
        pred = cdir / f"paired_{rung}_predictions.bin"
        recomputed, _po = failure_mask_from_predictions(cdir / "shots_bposd.bin", pred)
        masks[rung] = recomputed
        check(f"ladder.paired_failures_recomputed[{rung}]",
              int(recomputed.sum()) == int(paired["failures"][rung]),
              f"{int(recomputed.sum())} of {shots}")
    names = list(masks)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            key = f"{a}|{b}"
            if key not in paired["discordance"]:
                continue
            fa, fb = masks[a], masks[b]
            n_ab, n_ba = int((fa & ~fb).sum()), int((fb & ~fa).sum())
            rec = paired["discordance"][key]
            ok = (n_ab == rec["a_fails_b_ok"] and n_ba == rec["b_fails_a_ok"]
                  and abs(mcnemar_exact(n_ab, n_ba) - rec["mcnemar_exact_p"]) < 1e-12)
            check(f"ladder.discordance_recomputed[{key}]", ok,
                  f"{n_ab}/{n_ba} mcnemar_p={rec['mcnemar_exact_p']:.4g}")
    return summary


def verify_mechanism(cdir: Path) -> dict:
    summary = json.loads((cdir / "summary.json").read_text())
    manifest = json.loads((cdir / "manifest.json").read_text())
    shots = int(summary["shots_shared"])
    check("mech.campaign_closed", (cdir / "inventory.json").exists(), cdir.name)
    check_frozen_hashes("mech", cdir)
    check("mech.clean_3e-3_artifact",
          manifest["dem_dimensions"]["observables"] == 12
          and "argrescale" in manifest["circuit_file"],
          f"{manifest['circuit_file']} obs={manifest['dem_dimensions']['observables']}")
    eq = manifest["beam_equivalence_evidence"]
    check("mech.equivalence_gate",
          len(eq) == 6 and all(item["mismatched_shots"] == 0 for item in eq),
          f"{len(eq)} reports, all zero-mismatch")
    sample = json.loads((cdir / "sample_inventory.json").read_text())["shared_stream"]
    check("mech.sample_hash",
          sha256_file(cdir / "shots_shared.bin") == sample["sha256"],
          sample["sha256"][:12])
    check("mech.sampling_seed",
          manifest["sampling_seed"] == sample["sampling_seed"],
          str(sample["sampling_seed"]))
    masks = {}
    for arm in ("beam8", "beam32", "beam64"):
        recomputed, _po = failure_mask_from_predictions(
            cdir / "shots_shared.bin", cdir / f"{arm}_predictions.bin")
        masks[arm] = recomputed
        check(f"mech.failures_recomputed[{arm}]",
              int(recomputed.sum()) == int(summary["failures"][arm]),
              f"{int(recomputed.sum())} of {shots}")
    shard_data = [json.loads(p.read_text()) for p in sorted(cdir.glob("bposd_shard_*.json"))]
    bp = np.concatenate([read_failmask(cdir / s["failure_mask"], int(s["shots"]))
                         for s in sorted(shard_data, key=lambda s: int(s["shard"]))])
    masks["bp30+osd"] = bp
    check("mech.failures_recomputed[bp30+osd]",
          int(bp.sum()) == int(summary["failures"]["bp30+osd"]), f"{int(bp.sum())} of {shots}")
    f8, f32, f64 = (int(masks[k].sum()) for k in ("beam8", "beam32", "beam64"))
    check("mech.M1_ordering", f8 >= f32 >= f64, f"{f8} >= {f32} >= {f64}")
    check("mech.M3_ladder_sanity", f64 <= f32, f"{f64} <= {f32}")
    check("mech.clauses_match_summary",
          summary["clauses"]["M1_expansion_monotonicity"]["ordering_ok"] == (f8 >= f32 >= f64)
          and summary["clauses"]["M3_ladder_sanity"]["pass"] == (f64 <= f32),
          summary["mechanism_verdict"])
    return summary


def verify_gb7(cdir: Path) -> dict:
    summary = json.loads((cdir / "summary.json").read_text())
    manifest = json.loads((cdir / "manifest.json").read_text())
    sample = json.loads((cdir / "sample_inventory.json").read_text())
    shots = int(summary["shots_shared"])
    check("gb7.campaign_closed", (cdir / "inventory.json").exists(), cdir.name)
    check_frozen_hashes("gb7", cdir)
    check("gb7.prefix_campaign_pinned",
          manifest["prefix_campaign"]["run_id"] == PREFIX_CAMPAIGN.name,
          manifest["prefix_campaign"]["run_id"])
    check("gb7.binary_hash",
          manifest["beam_binary_sha256"] == sha256_file(BEAM_BINARY),
          manifest["beam_binary_sha256"][:12])
    eq = manifest["beam_equivalence_evidence"]
    check("gb7.equivalence_gate",
          len(eq) == 6 and all(item["mismatched_shots"] == 0 for item in eq),
          f"{len(eq)} reports, all zero-mismatch")
    shot_path = cdir / "shots_bposd.bin"
    check("gb7.sample_hash",
          sha256_file(shot_path) == sample["shared_stream"]["sha256"],
          sample["shared_stream"]["sha256"][:12])
    shot_prefix = compare_shot_prefix(
        shot_path, PREFIX_CAMPAIGN / "shots_bposd.bin")
    check("gb7.shot_prefix_continuity",
          all(record["byte_match"] for record in shot_prefix.values()),
          "detectors + observables")

    masks: dict[str, np.ndarray] = {}
    for rung in RUNGS:
        prediction = cdir / f"{rung}_predictions.bin"
        arm = summary["arms"][rung]
        check(f"gb7.predictions_hash[{rung}]",
              sha256_file(prediction) == arm["decode"]["predictions_sha256"],
              arm["decode"]["predictions_sha256"][:12])
        prefix = compare_prediction_prefix(
            prediction,
            PREFIX_CAMPAIGN / f"paired_{rung}_predictions.bin",
            int(manifest["dem_dimensions"]["observables"]),
        )
        check(f"gb7.prediction_prefix[{rung}]",
              prefix["byte_match"] is True, prefix["sha256"][:12])
        recomputed, _per_observable = failure_mask_from_predictions(
            shot_path, prediction)
        masks[rung] = recomputed
        check(f"gb7.failures_recomputed[{rung}]",
              int(recomputed.sum()) == int(arm["failures"]),
              f"{int(recomputed.sum())} of {shots}")
        frozen_mask = read_failmask(
            PREFIX_CAMPAIGN / f"paired_{rung}_failmask.bin", PREFIX_SHOTS)
        check(f"gb7.failure_mask_prefix[{rung}]",
              np.array_equal(recomputed[:PREFIX_SHOTS], frozen_mask),
              f"{int(recomputed[:PREFIX_SHOTS].sum())} prefix failures")

    for rung in ("beam32", "beam64"):
        redone = rung_decision(
            rung, paired_cells(masks["beam8"], masks[rung]), shots)
        stored = summary["decisions"][rung]
        same = (
            redone["outcome"] == stored["outcome"]
            and redone["paired_cells"] == stored["paired_cells"]
            and redone["paired_ratio_rung_over_beam8"]["ci95"]
            == stored["paired_ratio_rung_over_beam8"]["ci95"]
            and redone["mcnemar"]["exact_two_sided_p"]
            == stored["mcnemar"]["exact_two_sided_p"]
        )
        check(f"gb7.decision_reproducible[{rung}]", same,
              f"{redone['outcome']} ratio="
              f"{redone['paired_ratio_rung_over_beam8']['ratio']} "
              f"CI={redone['paired_ratio_rung_over_beam8']['ci95']}")
    return summary


def verify_gb8(gb7_summary: dict, gb7_dir: Path) -> dict:
    """Revision GB8: the re-target artifact must be a pure function of the
    frozen GB7 summary, and the corrected beam64 target must sit inside the
    frozen paired interval (the documented CONSISTENT verdict)."""
    check("gb8.artifact_present", GB8_ARTIFACT.is_file(), GB8_ARTIFACT.name)
    stored = json.loads(GB8_ARTIFACT.read_text())
    redone = retarget_all(gb7_summary)
    check("gb8.artifact_reproducible",
          all(stored[k] == redone[k] for k in ("rungs", "paper", "not_tested",
                                              "gb7_campaign")),
          f"from {gb7_dir.name}")
    check("gb8.source_pinned",
          stored["source"]["campaign"] == gb7_dir.name
          and stored["source"]["summary_sha256"] == sha256_file(gb7_dir / "summary.json"),
          stored["source"]["summary_sha256"][:12])
    b64 = redone["rungs"]["beam64"]
    check("gb8.beam64_paper_row_is_num_results_1",
          b64["paper_configuration"] == "beam64_640iters"
          and b64["paper_parameters"]["num_results"] == 1
          and GB8_PAPER["section3_factors_vs_bposd"]["beam64_640iters"] == 7.0,
          "7.0x, not the 17x of beam64_32res_640iters")
    corrected = b64["corrected_target"]
    lo, hi = b64["frozen_gb7"]["ci95"]
    check("gb8.beam64_corrected_target_inside_ci95",
          lo <= corrected["published_reciprocal"] <= hi
          and corrected["outcome"] == "NO_GAP_WIDTH_EFFECT_CONFIRMED_AND_INTERVAL_EXCLUDES_ONE",
          f"1/7.0={corrected['published_reciprocal']:.4f} in [{lo:.4f}, {hi:.4f}]")
    check("gb8.as_run_target_replay_matches_frozen",
          all(redone["rungs"][r]["gb7_target_as_run"]["outcome"]
              == gb7_summary["decisions"][r]["outcome"] for r in ("beam32", "beam64")),
          "beam32 NO_GAP, beam64 GAP_CONFIRMED")
    return redone


def main() -> None:
    ladder_dir = Path(sys.argv[1]).resolve()
    ladder = verify_ladder(ladder_dir)
    mech = verify_mechanism(Path(sys.argv[2]).resolve()) if len(sys.argv) > 2 else None
    gb7 = verify_gb7(Path(sys.argv[3]).resolve()) if len(sys.argv) > 3 else None
    gb8 = verify_gb8(gb7, Path(sys.argv[3]).resolve()) if gb7 is not None else None
    ps = (TARGET / "pre_statement.md").read_text()
    for rev in ("Revision GB5 ", "Revision GB5a ", "Revision GB6 ", "Revision GB7 ",
                "Revision GB8 "):
        check(f"pre_statement.contains[{rev.strip()}]", rev in ps)

    print()
    if failures:
        print(f"QLDPC_GATEB_ACCEPTANCE_FAIL {len(failures)} checks failed: {failures}")
        raise SystemExit(1)
    payload = {
        "ladder": ladder["ladder_verdict"],
        "rung_ratios": {r: {"ratio": v.get("beam_over_bposd_ratio"), "ci": v["ratio_ci95"],
                            "outcome": v["outcome"]}
                        for r, v in ladder["rung_verdicts"].items()},
        "paired_failures": ladder["paired_instrument"]["failures"],
    }
    if mech is not None:
        payload["mechanism"] = {"failures": mech["failures"],
                                "verdict": mech["mechanism_verdict"],
                                "clauses": {k: v["pass"] for k, v in mech["clauses"].items()}}
    if gb7 is not None:
        payload["gb7"] = {
            "failures": {rung: arm["failures"] for rung, arm in gb7["arms"].items()},
            "outcomes": {rung: decision["outcome"]
                         for rung, decision in gb7["decisions"].items()},
        }
    if gb8 is not None:
        payload["gb8"] = {
            rung: {"published_factor": r["corrected_target"]["published_factor"],
                   "outcome": r["corrected_target"]["outcome"]}
            for rung, r in gb8["rungs"].items()
        }
    print("QLDPC_GATEB_ACCEPTANCE_PASS " + json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    main()
