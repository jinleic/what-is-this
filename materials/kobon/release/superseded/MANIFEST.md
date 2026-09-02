# Superseded release archives — `math/kobon/release/superseded/`

None of these files is deleted. Each was relocated here on 2026-09-01 from
`math/kobon/release/` after the canonical release was established from
documentary evidence. Restore any file with the exact command recorded in its
row (run from the workspace root `/Users/jinleic/jinleic-workspace`).

## Canonical ruling

- **Canonical release:** `math/kobon/release/kobon-2026-08.zip`
  (19,286,189 bytes, SHA-256 `51e3789d1bb7b571302e6ab990f5a31d0b21790a316484cb83fdd2a984fca70d`,
  matching sidecar `math/kobon/release/kobon-2026-08.zip.sha256`), together
  with its unpacked mirror `math/kobon/release/kobon-2026-08/`.
- Evidence: `math/kobon/docs/INDEX.md` §"Canonical release snapshot" ("The
  canonical archive is `math/kobon/release/kobon-2026-08.zip` … earlier
  `*-obsolete-*` and `-v1` archives are superseded packaging history, not
  alternate authorities"); `math/kobon/docs/ARTIFACTS.md` rows 92–96;
  `math/PROGRESS.md` entries of 2026-08-21 ("canonical archive
  `math/kobon/release/kobon-2026-08.zip` … 44/44 manifest verified"; "release
  v2 zip rebuilt (45/45 manifest verified)"); the v8 package README itself
  ("Companion canonical zip: `../kobon-2026-08.zip`").

## Version-chain finding

The v1→v8 series is a chain of packaging iterations of the same consolidated
snapshot (inner member path `kobon-2026-08/...`, growing 53 → 104 members,
v8's manifest has 104 entries), each superseding the previous as curation
progressed; the chain terminates in v8, after which the snapshot was repacked
as the final top-level canonical zip `kobon-2026-08.zip` (134 members, flat
layout). The `obsolete-structure*` trio are the earliest packaging drafts.
No versioned zip introduces results absent from the canonical zip: the
canonical zip's packaged `MANIFEST.sha256` is byte-identical (SHA-256
`8a464ea78bc7af19d9a4754b8914495c98a8730b511028273129ef9bab38c03d`) to the
live mirror's, which passes `sha256sum -c` with 133/133 OK.

## Integrity

All 12 files were hash-verified immediately before and after the move;
SHA-256 and byte size are unchanged for every file. Byte-identical duplicate
groups: **none** — all 13 release zips have distinct SHA-256 digests.

## Roster (12 files, 5,766,585 bytes total)

| original path (under `math/kobon/release/`) | bytes | mtime | SHA-256 | why superseded | restore command |
|---|---|---|---|---|---|
| `kobon-2026-08-obsolete-structure.zip` | 248,543 | 2026-08-21 19:25:07 UTC | `052485125b30f3cc31d2ccdd57f9a65d977936ae2aadb5a13f7d61018b3ad3a3` | earliest packaging draft of the consolidated release; self-labelled obsolete | `mv math/kobon/release/superseded/kobon-2026-08-obsolete-structure.zip math/kobon/release/` |
| `kobon-2026-08-obsolete-structure2.zip` | 248,571 | 2026-08-21 19:26:48 UTC | `bcd90fbcdffd32155d3df701d5568b03e1dd9e2c1ac912ebd0836935d49da2df` | second packaging draft; self-labelled obsolete | `mv math/kobon/release/superseded/kobon-2026-08-obsolete-structure2.zip math/kobon/release/` |
| `kobon-2026-08-obsolete-structure3.zip` | 248,513 | 2026-08-21 19:27:58 UTC | `f29ed2aa4f8a6cbae19a7ebf2c1cd60a60043265b0a132062ce036450d2cb822` | third packaging draft; self-labelled obsolete | `mv math/kobon/release/superseded/kobon-2026-08-obsolete-structure3.zip math/kobon/release/` |
| `kobon-2026-08-v1.zip` | 248,593 | 2026-08-21 19:30:16 UTC | `a4dfa51a5c728a49f21fe7c95e8572eaad08b0bca3dd235b13d3d5657c292f63` | version 1 of the wrapped-structure series; superseded by v2+ | `mv math/kobon/release/superseded/kobon-2026-08-v1.zip math/kobon/release/` |
| `kobon-2026-08-v2.zip` | 254,938 | 2026-08-21 20:58:48 UTC | `698b6c8136c563dbb63d856d76d02a2c851f3304db562ed52daff1821991270a` | version 2 (45/45 manifest); superseded by v3+ | `mv math/kobon/release/superseded/kobon-2026-08-v2.zip math/kobon/release/` |
| `kobon-2026-08-v3.zip` | 594,945 | 2026-08-22 02:25:06 UTC | `a0982f2455de2c8b29184fff4f93ad9820b92af272767ee89b965f4f13884548` | version 3 of the wrapped-structure series; superseded by v4+ | `mv math/kobon/release/superseded/kobon-2026-08-v3.zip math/kobon/release/` |
| `kobon-2026-08-v4.zip` | 616,111 | 2026-08-22 03:06:40 UTC | `75fca6e220383a97c3f2b1d5ce598d62b18ae27839541a8c1f6014bb1406ad41` | version 4 of the wrapped-structure series; superseded by v5+ | `mv math/kobon/release/superseded/kobon-2026-08-v4.zip math/kobon/release/` |
| `kobon-2026-08-v5.zip` | 618,892 | 2026-08-22 03:10:05 UTC | `5a345115741b32525aede177e3fcd00fbccc36ad3e2d04207d39006be6c657ae` | version 5 of the wrapped-structure series; superseded by v5b+ | `mv math/kobon/release/superseded/kobon-2026-08-v5.zip math/kobon/release/` |
| `kobon-2026-08-v5b.zip` | 619,154 | 2026-08-22 03:11:13 UTC | `9b0d490006bb137dbaadf32af565d584e2d5cdb495333bbb021b0f4fd5c15cab` | v5 correction; superseded by v6+ | `mv math/kobon/release/superseded/kobon-2026-08-v5b.zip math/kobon/release/` |
| `kobon-2026-08-v6.zip` | 661,330 | 2026-08-22 05:48:34 UTC | `1926b6d3fc1328611d04077a7b2fa331f7b50a18d7fcad0c3e026910eb5a7fa2` | version 6 of the wrapped-structure series; superseded by v7+ | `mv math/kobon/release/superseded/kobon-2026-08-v6.zip math/kobon/release/` |
| `kobon-2026-08-v7.zip` | 679,926 | 2026-08-22 06:23:00 UTC | `282621b5fe47d31155270c4cfecd734540d6651a8b392fe3506a2ac4b21adfbd` | version 7 of the wrapped-structure series; superseded by v8 | `mv math/kobon/release/superseded/kobon-2026-08-v7.zip math/kobon/release/` |
| `kobon-2026-08-v8.zip` | 727,069 | 2026-08-22 20:51:43 UTC | `85d2f01f22e03be6f97d3c62a0ad040ac8fdc9f72a312b5cc2c8893561b0e8d4` | last wrapped-structure iteration; its own README defers to `../kobon-2026-08.zip` as canonical | `mv math/kobon/release/superseded/kobon-2026-08-v8.zip math/kobon/release/` |

## Provenance of this ruling

Supersession is documentary, not inferential: `INDEX.md` names the canonical
archive explicitly; `ARTIFACTS.md` labels v1 and `obsolete-structure*` rows
"referenced only (superseded archive)"; `PROGRESS.md` narrates the packaging
walk (44/44 → 45/45 manifests) on 2026-08-21; the v8 package README
defers to the canonical zip. The 12 moved files are still referenced by name
inside `math/kobon/docs/ARTIFACTS.md` snapshot tables (point-in-time
inventory, no path dependency), and nowhere else in the workspace.
