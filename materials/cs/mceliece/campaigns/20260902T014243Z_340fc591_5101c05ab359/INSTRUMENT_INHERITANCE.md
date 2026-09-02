# Instrument inheritance record — t=48 fresh evidential run (recorded at INIT, before ANY compute)

Fresh post-t=96 evidential campaign (`20260902T014243Z_340fc591_5101c05ab359`),
init'd 2026-09-02T01:42:43Z from the COMMITTED prereg `7c05758`
(`cs/mceliece/pre_statement_t48.md`, sha256
0ab39d4bc10f70e80a1a506ab06947bda27ad0927bb7e484688ab57e7ba9c867 —
confirmed by init itself). Registered cell (12, 3488, 48, 16384),
k=2912, D=3391, 2t+3=99.

Instrument inheritance, verified in-place before any t=48 compute:

- `code/{instance,fastfield,gfield,census}.py` BYTE-IDENTICAL to the
  current `src/` copies (`cmp -s` exit 0, all four).
- BYTE-IDENTICAL to the certified `M12ALPHA_DIVONLY` frozen snapshot,
  all four.
- BYTE-IDENTICAL to the certified t=96 run's copies
  (`../20260901T100025Z_bd55fdfd_8fbba9be5ab3/code/`), all four.
- (Historic third reference, the 14-09Z ladder snapshot: differs only
  by the inert `forced_G` ctor param, as documented by both prior
  campaigns — same file set, same position.)

Stage scripts in THIS run `verdict_m12_t48_stage.py`,
`anchors_stage.py`, `plants_m12_t48_stage.py` are mechanical
re-instantiations of the certified t=96 stages: rename-equivalence
diffs recorded (only doc-title case, the run id, the cell tuple, and
t96->t48 variable-name spellings differ; logic untouched). The plants
stage's P-4 probe now points at (12,3488,48,16384) per pre-statement
section 4.

Both prior m=12 runs (t=96 FROZEN-CERTIFIED, its REHEARSAL
predecessor) are untouched; this campaign neither re-runs nor extends
them.

— McelieceM12T96, 2026-09-02T01:4xZ, at init of the t=48 run, before any compute.
