# RLIMIT diagnosis + run-log record (2026-09-02, before relaunch)

Main diagnosis adopted: exit 152 = 128+24 = SIGXCPU (CPU-time rlimit),
NOT OOM. Evidence: first verdict run died after ~10.7 min wall; my
launch context shows `resource.getrlimit(resource.RLIMIT_CPU)` =
(soft=7200 s, hard=RLIM_INFINITY). 7200 s soft ≈ kill after ~2 h CPU —
matches the second nohup run dying silently around the 2 h mark (the
log stayed empty because the stage writes only at build completion;
outer `timeout` accounting made the first bg_6 death look like ~10 min
wall but ~600-7200 s CPU class). The t=64/t=96 verdict stages ran in
the same kernel-spawned contexts with the same rlimit and never
exceeded their (shorter) CPU needs, which is why this never surfaced.
 swap pressure note retracted as the cause.

Mitigation: raise the soft limit to the hard limit
(resource.setrlimit(RLIMIT_CPU, (hard, hard))) at the top of the
relaunched process via a bash `python3 -c` wrapper — the stage script
itself stays UNCHANGED (no amendment needed; Main directive). nice
value disclosed: 10 (box load is a peer's beam8_cpp at nice 10; per
Main, avoid 15; instrument loop discipline unchanged — 'nice -n 10'
was the pre-registered instrumentation anyway).

Recorded values: RLIMIT_CPU in bash-tool launch context:
soft=7200 hard=9223372036854775807 (RLIM_INFINITY).

— McelieceM12T96
