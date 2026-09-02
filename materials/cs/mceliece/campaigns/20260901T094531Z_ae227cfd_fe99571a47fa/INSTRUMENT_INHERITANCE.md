# Instrument inheritance record — byte-diff vs frozen snapshots (recorded at INIT, before any compute)

The inherited instrument is the certified `M12ALPHA_DIVONLY` code path,
copied byte-for-byte into `code/` at campaign init
(2026-09-01T09:45:31Z), before any t=96 compute:

- `instance.py`, `fastfield.py`, `gfield.py`, `census.py` in `code/`
  are BYTE-IDENTICAL to the current `src/` copies
  (`cmp -s` exit 0 for all four; recorded in the shell transcript).
- `src/` itself is byte-identical to the `M12ALPHA_DIVONLY` frozen
  snapshot `campaigns/2026-09-01T03-28-10Z_M12ALPHA_DIVONLY/code/`
  (all four files, `cmp -s` exit 0).
- vs the frozen `2026-08-30T14-09Z_57200ADD/code/` snapshot:
  `fastfield.py`, `gfield.py`, `census.py` BYTE-IDENTICAL;
  `instance.py` differs ONLY in lines 84–91 — the inert `forced_G`
  constructor parameter (an unused-input artifact of Gate-B tooling,
  disclosed by the t=64 campaign; present in the current file, absent
  in the 14-09Z snapshot). Same diff the t=64 campaign recorded.

Conclusion: THIS campaign's t=96 verdicts run on exactly the instrument
that produced the certified m≤11 ladder (14-09Z), the Gate-B rate
verifier, and the certified (12,3488,64,16384) row. NO second
convention is created. The delta/alpha stage script
(`verdict_m12_t96_stage.py`) is a mechanical re-instantiation of
`M12ALPHA_DIVONLY/verdict_m12_stage.py` at the new registered cell —
diffed against it at freeze time; the anchors and plants scripts
likewise inherit the t=64 campaign's `camp_lib.py` /
`plants_stage.py` logic with the cell tuple changed.

— McelieceM12T96, 2026-09-01, at init, before compute.
