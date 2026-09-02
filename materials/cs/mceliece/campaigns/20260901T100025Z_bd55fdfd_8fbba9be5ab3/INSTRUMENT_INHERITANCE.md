# Instrument inheritance record — fresh evidential run (recorded at INIT, before ANY compute)

This is the fresh post-rehearsal evidential campaign
(`20260901T100025Z_bd55fdfd_8fbba9be5ab3`), init'd at 2026-09-01T10:00:25Z
from the COMMITTED prereg `9189ff7`
(`cs/mceliece/pre_statement_t96.md`, sha256
9d1e28d56c9e9af777c0bae062b375607e2a7857be599c4e37071f4e800d3047 —
confirmed byte-identical by init itself).

Instrument inheritance, verified in-place before any t=96 compute:

- `code/{instance,fastfield,gfield,census}.py` are BYTE-IDENTICAL to
  the current `src/` copies (`cmp -s` exit 0, all four).
- BYTE-IDENTICAL to the certified `M12ALPHA_DIVONLY` frozen snapshot
  (`campaigns/2026-09-01T03-28-10Z_M12ALPHA_DIVONLY/code/`), all four.
- vs the frozen `2026-08-30T14-09Z_57200ADD/code/` ladder snapshot:
  `fastfield.py`, `gfield.py`, `census.py` BYTE-IDENTICAL;
  `instance.py` differs ONLY at lines 84–91 — the inert `forced_G`
  constructor parameter (unused-input artifact of Gate-B tooling,
  disclosed by the t=64 campaign). Identical diff position/content to
  what both prior campaigns recorded.

Stage scripts for THIS run are created after init (dated by mtime),
mechanically re-instantiating the M12ALPHA_DIVONLY stages at the
registered t=96 cell; they will be frozen with the run. No second
convention is created.

Preceding rehearsal run (`20260901T094531Z_ae227cfd_fe99571a47fa`,
closed REHEARSAL) and ALL of its outputs are non-evidential and are
NOT reused; its provenance: `REHEARSAL_PROVENANCE.md` in that run dir.

— McelieceM12T96, 2026-09-01T10:01Z, at init of the fresh run, before any compute.
