# Session → project map

omp session IDs in `~/.omp/agent/sessions/-jinleic-workspace/<timestamp>_<id>.jsonl`
(head record: title + cwd; first user message names the problem).

| Session | Started | Problem |
|---|---|---|
| `019ff0a0-6557-7000-8b8c-9bdbccb719e4` | 2026-08-11 | union-closed (uc); Millennium/NS prelude; ccf subagents |
| `019ff186-f657-7000-85d0-6d52723bfb5c` | 2026-08-11 | 3D Ising (ising3d) |
| `019ff1a6-af1a-7000-b0e9-223c45ae0a3c` | 2026-08-11 | quantum LDPC / QEC |
| `019ff95f-d695-7000-b3d0-36dd86537081` | 2026-08-13 | Kobon triangles |
| `019ff957-1ed4-7000-8941-f92f3002c676` | 2026-08-13 | Hilbert's tenth over Q (h10q) |
| `019ffb23-9512-7000-8c79-32667475fc6c` | 2026-08-13 | R(5,5) + zeta(5) triage (r55, zeta5) |

Note: `ns` and `ccf` have no dedicated root session (ns was the opening scan of
`019ff0a0`; ccf ran as subagents inside it). Restored 2026-08-15 from the
migration bundle; audit: `SessionIndex` — `019ff0a0` is a union-closed cluster.
