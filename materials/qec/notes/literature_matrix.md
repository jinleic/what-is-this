# Literature matrix

Date checked: 2026-08-11

| Paper | What it establishes | What it leaves open for this project |
|---|---|---|
| Bravyi et al., arXiv:2308.07915 v2 | CSS bivariate-bicycle memory with 0.8% circuit-level threshold, depth-7 syndrome cycle, `n` ancillas, degree-6 connectivity, and the reported 12-logical-qubit memory result. | Circuit-level behavior of the non-CSS PBB perturbations. |
| Cruz-Benito et al., arXiv:2606.02418 v1 | Five campaigns; approximately 1650 iterations; 2e5 candidates; 465 codes (97 CSS, 368 PBB); selected CSS/PBB codes; term PBB; approximately 140h / approximately $400. | RESOLVED: the commutation condition is symmetry of `A C^T + B D^T`, stated correctly in Sec. III.2 and independently re-derived by us -- there is no paper/code discrepancy; our own paraphrase was at fault (`notes/failed_routes.md` FR-006). STILL OPEN: exact distances across the whole catalogue, and circuit-level decoding performance beyond the `[[144,12,12]]` family measured here. |
| deMarti iOlius et al., arXiv:2409.01440 v3 | BP+OTF decoder for quantum LDPC codes under circuit-level noise. | PBB-specific circuit-level comparison. |
| Voss et al., arXiv:2406.19151 v4 | Multivariate Bicycle Codes. | The PBB circuit-level and catalogue questions studied here. |
| Postema and Kokkelmans, arXiv:2502.17052 v4 (rev 2026-04-30) | Existence and characterisation of bivariate bicycle codes. | The PBB perturbation and its circuit-level behavior. |
| Yoder et al., arXiv:2506.03094 v1 | Tour de gross modular quantum computer based on bivariate bicycle codes. | Non-CSS PBB circuit-level results. |
| Mayer et al., arXiv:2509.13678 v2 | Rare-event simulation of quantum error-correction circuits. | PBB-specific rare-event estimates. |
| Chen et al., arXiv:2501.14380 v3 | Verification of fault tolerance for quantum error-correction codes. | Verification of the PBB circuits and decoders used here. |
| Palsberg et al., arXiv:2601.20247 v1 | Computer science challenges in quantum computing, including early fault tolerance and beyond. | The specific PBB mathematical and circuit-level questions here. |
| Webster, Jacob, and Higgott, arXiv:2603.22532 v1 | Exact distance methods versus heuristic distance methods, including SAT/CLISAT, MIP-SCIP/Gurobi, BP-OSD, and Stim methods. | Exact distances for the remaining catalogue codes and PBB circuit-level results. |
| Chengyu et al., arXiv:2601.18562 v1 | Bayesian optimization for quantum error-correcting code discovery; repository is identified as `https://github.com/Chengyuyihua/BO_for_QECcodes`. | Relationship to the PBB catalogue and circuit-level performance. |
| QDistRnd, arXiv:2308.15140 | GAP package for computing quantum-code distance. | Exact distances of the remaining PBB catalogue codes by independent methods. |
| Khesin and Lu, arXiv:2603.05496 (March 2026) | SECONDHAND / needs direct verification: Mirror codes are a non-CSS generalization of abelian two-block group-algebra codes; reported best weight-6 non-CSS examples are `[[60,4,10]]`, `[[36,6,6]]`, `[[48,8,6]]`, `[[85,8,9]]`. | Direct verification and comparison with PBB results. |
