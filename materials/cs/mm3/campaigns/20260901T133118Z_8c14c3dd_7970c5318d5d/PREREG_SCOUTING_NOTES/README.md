# Pre-registration scouting record (NON-EVIDENTIAL)

These scripts ran BEFORE the prereg commit 9a90db0 to size the budget and
pin the universes. None ran under campaign.py; no CNF/kissat instance was
run pre-registration; they are disclosed in pre_statement.md §8 and are
NOT evidence for any verdict. Every claimed campaign result is re-derived
in-run after init under the dual-instrument protocol from scratch.

Recorded scouting observations (instrument 1 = frozen subset-DFS only):
- mws59 floor re-runs: d=(13,12,14), floor impossible, states (132,81,672).
- U-side aux-1 full DFS scan at T=14: 1 of 389 feasible — a=(0,0,0,-1,1,1,-1,0,1).
- V-side aux-1 full DFS scan at T=13: 0 of 366 feasible.
- W-side aux-1 full DFS scan at T=15: 8 of 421 feasible.
- V-side pair-DFS timed on ~130 sampled pairs: 0 feasible in sample;
  48-65 ms/pair measured.
- Universe hashes: tau(U)=13 f464b756…, tau(V)=12 824c0f0e…, tau(W)=14
  2d860371…; AU(U)=389 b0d0584e…, AU(V)=366 9b060604…, AU(W)=421 5aaa81dc….
- Timing scouts: aux-1 DFS 35-65 ms/instance; pair DFS ~52 ms/pair;
  DFS peak RSS 17.2 MiB.
