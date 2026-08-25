# Artifact inventory

Inventory base snapshot: **2026-08-22T01:55:20Z**. Paths are workspace-relative. The large tables below preserve that point-in-time campaign inventory. The 2026-08-23 square-penalty frontier delta is listed separately below; for the current release contents and hashes, `math/kobon/release/kobon-2026-08/MANIFEST.sha256` is authoritative. Third-party virtual environments, solver binaries/source trees, `__pycache__`, transient process markers, and uncited superseded exploratory CNFs are dependencies or scratch state rather than campaign artifacts and are excluded.

Files strictly smaller than **50,000,000 bytes** have a SHA-256 computed from disk. Larger files record `hash omitted (large)`, as requested. Hashes of in-flight logs are point-in-time values and are expected to change. Zip inclusion is determined against `math/kobon/release/kobon-2026-08/MANIFEST.sha256`; “included byte-identically” may name the path used inside the zip. `ARTIFACTS.md` cannot contain a stable hash of itself, so it is the sole self-referential omission.

**Inventory totals:** 427 files; 302 SHA-256 values computed; 125 large-file hashes omitted.

## 2026-08-23/24 frontier deltas

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/square_penalty_counterexample.py` | 14,230 | `bd695af8fd3b6f441fd7bce8780159f4751b2f473ea8b426b209500cc6d7684d` | included as `scripts/square_penalty_counterexample.py` |
| `scratch/kobon/square_penalty_counterexample.json` | 7,237 | `6555f65b47dfc3f5d08521fa9a99e97f4193f3f5309a70d867a8d3efaecd95a9` | included as `verification/square_penalty_counterexample.json` |
| `scratch/kobon/square_penalty_sat.py` | 6,759 | `43a676c756d67c939417fc582d4ca065a91896d4ad4e03b43c139a6421aee0e9` | included as `scripts/square_penalty_sat.py` |
| `scratch/kobon/c23_n7_violation.cnf` | 2,533,909 | `075ca17768ad70e4ff69673543c85d629a69737656b36211c949b8e0a32676d9` | included as `certificates/c23_n7_violation.cnf` |
| `scratch/kobon/c23_n7_violation.drat` | 34,133,280 | `b17bb8d6610242ee3149b141d24447b2261ab1c75ac2f4e92c21913ad53f4c5a` | included as `certificates/c23_n7_violation.drat` |
| `scratch/kobon/c23_n7_violation.dratcheck.log` | 442 | checker transcript ending `s VERIFIED` | included as `certificates/c23_n7_violation.dratcheck.log` |
| `scratch/kobon/c23_n7_certificate.json` | 2,085 | `214c5fca509c4de1a310385e8d63cc0bf0e47ef25dfd57446a548f8ad6912b9b` | included as `certificates/c23_n7_certificate.json` |
| `math/kobon/engine.py` | 65,426 | `d1ff89892e0bc5130c701365066c0cd43dbdc3d4d1d96e836502115617ff6351` | included as `scripts/engine.py` |
| `math/kobon/test_engine.py` | 32,738 | `ebcc5feaec5d2f582ba28504f2590f935f1eacb3ba2f9d8acc3696090a166885` | referenced only; focused suite 36/36 |
| `scratch/kobon/gap_faces.py` | 8,859 | `4b2f9668993b7be5cf56135b53a430990f37c19aacf68ab79df7904c7cf1c412` | included as `scripts/gap_faces.py` |
| `scratch/kobon/sector_bound_audit.py` | 14,458 | `d4978d7b8a7893245a8d4dbbd0c332ac30116b67259498b07694c75b5524c6a9` | included as `verification/sector_bound_audit.py` |
| `scratch/kobon/sector_bound_audit.json` | 1,280 | `b2bd2fc173ae659a6cd1f16374903673a66094a4858c4a011fb491330bc56b04` | included as `verification/sector_bound_audit.json` |
| `scratch/kobon/pappus_relaxation_probe.py` | 3,601 | `975af524f5aad770d02dd0e20f22200f8054c9755d12821b8315f2bdd21d0880` | included as `verification/pappus_relaxation_probe.py` |
| `scratch/kobon/pappus_relaxation_probe.json` | 659 | `ba4a108a77f9f527fbd576c64d173d8f715aa5735ba6d68584dfbf5968c4550c` | included as `verification/pappus_relaxation_probe.json` |
| `scratch/kobon/n12_gap_sector_deletion_t39.cnf` | 16,145,451 | `5435cb9958710878b1a68bc112e353a46d9a031a723e1afaa31191a9a9546909` | referenced only; live discovery instance, no verdict |
| `scratch/kobon/n12_gap_sector_deletion_t39.metadata.json` | 3,727 | `47846eac53412110d8bf1d9053476e0e676378ab0a6807f8659fd570e4757d0b` | included as `discoveries/n12_gap_sector_deletion_t39.metadata.json` |
| `scratch/kobon/ladder_n10_t26_faces.cnf` | 22,125,971 | `d3eddc9c5ee996dc3946eeac1fdd3289daef8ffe0802199ca0d7a3791141dfa2` | included as `certificates/ladder_n10_t26_faces.cnf` |
| `scratch/kobon/ladder_n10_t26_faces.drat` | 14,050,806,440 | `bdbab47463806ee96163d03ec36a052e5818e2125202272182ee4cbe94629007` | referenced only; complete checked proof |
| `scratch/kobon/ladder_n10_t26_faces.dratcheck.log` | 452 | `d4efb6701373c33d259b8e08abe66b220c763b99e58d011578c915b380b22313` | included as `certificates/ladder_n10_t26_faces.dratcheck.log` |
| `scratch/kobon/ladder_n10_t26_faces_certificate.json` | 3,066 | `50085b4607ab70bc8e124472f77db1bb981a76266c9be18c1eb9de46a07795f9` | included as `certificates/ladder_n10_t26_faces_certificate.json` |
| `scratch/kobon/n12/n12_face_cells_q1_t39.cnf` | 13,615,953 | `fa8460252010e62bf9b722d9c08562e8ae6c2ddf996059901e9bce23e7c27bc4` | referenced only; discovery UNSAT without proof |
| `scratch/kobon/n12/n12_face_cells_q1_t39.discovery.json` | 2,875 | `3c64a0ed92edad42a3452b449ae9f8950cc85c85d7ee193c975d419a9d1ded70` | included as `discoveries/n12_face_cells_q1_t39.discovery.json` |
| `scratch/kobon/n12_gap_shared_ray_deletion_t39.cnf` | 16,474,653 | `0aa9e81e1806d5ff7413669d7357454c9e08a496d2172d8b4a512b27115b7e7d` | referenced only; live opt-in discovery instance, no verdict |
| `scratch/kobon/shared_ray_bound_experiment.json` | 5,050 | `3d434eea77da033fcf562d1d3ac552a38821f440774f59c705615ff6fa562fdb` | included as `discoveries/shared_ray_bound_experiment.json` |
| `scratch/kobon/reified_face_audit.py` | 7,156 | `76e15a255072e3a1a14075b7cb39cc0cc99d2b54461bdd74049acd0be03688b7` | included as `verification/reified_face_audit.py` |
| `scratch/kobon/reified_face_audit.json` | 1,108 | `24dbaf97250b5750e415b2132aaeb9d4df12c1eb8ac9ec56f507fb79452026f8` | included as `verification/reified_face_audit.json` |
| `scratch/kobon/chirotope_gp_audit.py` | 7,713 | `67fa614dbdafb540d17572f097e1d2e135bd50941c908bd9396a5af02c039676` | included as `verification/chirotope_gp_audit.py` |
| `scratch/kobon/chirotope_gp_audit.json` | 662 | `65395eadc15ba3326df7f14753fbe7c0aa7c85a26c5702d47844ba1803dea9a2` | included as `verification/chirotope_gp_audit.json` |
| `scratch/kobon/frontier_endpoint_research.json` | 9,556 | `3cda3792387f49e54c56be828be9b30c1209d07fa28965ebb3dc658c1e05c4ed` | included as `discoveries/frontier_endpoint_research.json` |
| `scratch/kobon/n12_gap_endpoint_mi_k4_deletion_t39.cnf` | 21,685,650 | `6dca104aa120ee1f209aada8aca07a08bf3341f1335e6e99f4ae1f11bc2ba76b` | referenced only; live discovery instance, no verdict/proof |
| `scratch/kobon/n12_gap_endpoint_mi_k4_sub9_t39.cnf` | 39,495,970 | `34406004c932ce77d33f76f8b7d19c7b0f5a2e3e241e94651846269bf750f2df` | referenced only; live \(K(9)\)-hereditary discovery instance, no verdict/proof |
| `math/kobon/paper/kobon_broad_capacity.tex` | 67,966 | `e39b5db1eb181de2f6b31894c358218fd9a15a17bbd3030bd0d7dd70063d5e7d` | included as `papers/kobon_broad_capacity.tex` |
| `math/kobon/paper/kobon_broad_capacity.pdf` | 377,434 | `1f794df33d7adfc1e35146191238763f374b66725a9aa30d22b4bddb7198523b` | included as `papers/kobon_broad_capacity.pdf` |

## Canonical and working documents

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `math/PROGRESS.md` | 439,839 | `7f66e620d7d0ca300114a1f455d9e8a27dd6f44edcec33a7ed34ee5154589ef0` | included as `docs/PROGRESS.md` |
| `math/kobon/README.md` | 10,178 | `aab60eda88312f6ee0df305beb48f529ec1651022862eb2056707871efd13d73` | referenced only |
| `math/kobon/AUDIT.md` | 7,567 | `981628f74ca64a4d419f3826521e9a18d7ad4407ed2dd97b179b93216c0117a9` | referenced only |
| `math/kobon/report.md` | 24,616 | `d6b6e95ac27e1dac478f3dd3018ee7d929cae3862734379eedc0538a14877333` | included byte-identically as `papers/report.md` |
| `math/kobon/paper_capacity.md` | 70,546 | `b664cb49ca581bb5536df3065c9702843cfb2bcbd13dfd3f2f97fc4c01bd2da3` | included byte-identically as `papers/paper_capacity.md` |
| `math/kobon/paper_kobon_2026-08.md` | 9,973 | `63eeac4be91857da840487c2ef10bcd149c5bbb6d4f5735d8ea75e1c1ebeedef` | included byte-identically as `papers/paper_kobon_2026-08.md` |
| `math/kobon/docs/INDEX.md` | 19,733 | `091d8a580fc1d81a56e7f398d25bfe58ad6e2f753288fce2440c95b4f3f3e407` | included as `docs/INDEX.md` |
| `math/kobon/docs/methods.md` | 11,879 | `8dbd5ffba71a9702d1ef0396ff23afa250c4f2130a5bd7ea5482247591db221b` | included as `docs/methods.md` |
| `math/kobon/docs/experiments.md` | 15,240 | `4da9f824057c0f150fa5d691844c699e7a97b0763d4895aebfc0017b6732a752` | included as `docs/experiments.md` |
| `math/kobon/docs/reproduce.md` | 11,758 | `89e80faefce3b7f60d0ffcc404997ab49c7a1451c4b4d8a62c50f216c3c22eb4` | included as `docs/reproduce.md` |
| `math/kobon/verification/README.md` | 7,426 | `2894d5874fe5d410f4c6980e5b0b1fbd401de71555fb6c4ae45c58ffaa64f7a5` | included as `verification/README.md` |
| `scratch/kobon/litrefs.md` | 21,776 | `76a54b74a05d31d4db7641b767e3d8320893b6f36ff5561858890d187ca036d1` | included byte-identically as `docs/litrefs.md` |
| `math/kobon/docs/NEXT_BREAKTHROUGHS.md` | 7,072 | `2b9d0fe0dd52016d25a1ce99944264451b4248d56feb60e0f77b9ff8872f1ce3` | included as `docs/NEXT_BREAKTHROUGHS.md` |
| `scratch/kobon/n11_chain.md` | 12,187 | `bbec7665795d0aa699ed6ac5934322f8b4cae17c8c481650d3e1187e55233e51` | referenced only |
| `scratch/kobon/novelty_audit_capacity_theorem.md` | 13,153 | `a400184f72942a4ce598ccb7f351cfe016ce7e63bf9858d0e731b0aef9123731` | referenced only |
| `scratch/kobon/n14_face_status.md` | 7,859 | `cac2d77c34b11301f458c829111364a30ceb85d0aa5e7967a37ca981203c4eee` | referenced only |
| `scratch/kobon/n11/section_draft.md` | 3,887 | `79c7e658558b567f9915f8cd34138a994a5335c95b84fd18d312960c0ca6026a` | referenced only |
| `scratch/kobon/n11/equality_branches.md` | 9,879 | `b197b83ca48c2793ad3f4f89fc5039fd42bb26dc9b6662904358397a52843015` | referenced only |
| `scratch/kobon/n12/closeout_draft.md` | 9,905 | `db83b4b147473f354f05543e4869db1233e2b079e0e5e1d1b9c0e19ea497af33` | referenced only |
| `scratch/kobon/n12/path_arbiter_rec.md` | 12,805 | `1c64e2fe45929adff94d741540e9f329e55973f90efe42c41e6742267a4812ec` | referenced only |
| `scratch/kobon/n14/Kgen14_ge54_certificate.md` | 2,970 | `240fd7e18989719923fd7bdc0c328cf96052e86b9c6387110123c0f2cbf7dcfe` | included byte-identically as `proofs/Kgen14_ge54_certificate.md` |
| `scratch/kobon/n14/orbit_cube_plan.md` | 12,264 | `69a7b58b9de42a46097ad734251a903265c8b4bfb4006123b9ea897ca173aa1a` | included byte-identically as `docs/orbit_cube_plan.md` |
| `scratch/kobon/n14/escape_obstruction.md` | 16,862 | `def268d66bcd339ec3a61eaec66dd899106d4fb08c48bd90e2f079386bb7f424` | included byte-identically as `proofs/escape_obstruction.md` |
| `scratch/kobon/n14/concurrency_cubes.md` | 9,991 | `d1a50a7535f25996871689b45bb2392340cd8b504d06ba7247c7483ed3ddfa23` | referenced only |
| `scratch/kobon/n14/work/NOTES.md` | 5,289 | `0b14ab67d44ec98add40896a506f6468ae930e00455001e5a86e09ae0b803da0` | referenced only |
| `scratch/kobon/n14/work/line-order/README.md` | 8,785 | `1e51bb7ef771e87508542386d00aea87e00152248d564dce71dd7c238a8c16a2` | referenced only |
| `scratch/kobon/n14/work/line-order/CHANGELOG.md` | 1,546 | `a9f63743d6ec24ca290a037e98ffb41223cae8ff1a3b83dac256c2699c2d41de` | referenced only |
| `scratch/kobon/n14/work/line-order/TODO.md` | 1,269 | `9b819f1c1ecff7ff92daa35ff09477653afb153aa3751fbd066345621bc1d069` | referenced only |
| `scratch/kobon/n14/work/kobon-cnf/README.md` | 2,200 | `30142de1f162c95699d73f288c3dffdb0c6c7fca9734e6d742d72343c5f6581a` | referenced only |
| `scratch/kobon/n14/work/kobon-cnf/CHANGELOG.md` | 398 | `07ccf5f0a14daf365b3657ad45f1ad1d15391257ca29244c7ed01d30251830b5` | referenced only |
| `scratch/kobon/n14/work/kobon-cnf/TODO.md` | 183 | `c253d686d5e4598e88592041025cdf2bd3eb43b8b90c18b62c4b46393a4f22d6` | referenced only |
| `scratch/kobon/n18/obstruction_n18.md` | 18,686 | `5da2ddcfaedecab898c520ce6f08da55bde61dd510d75703660944b938ce20fa` | included byte-identically as `proofs/obstruction_n18.md` |
| `scratch/kobon/n18/batch2_launch.md` | 2,947 | `265635197a56125bfc81ffd5c96b9fcd7c8160f243c0111b553611c3252aa0f6` | included byte-identically as `docs/batch2_launch.md` |
| `scratch/kobon/n20/signature_algebra.md` | 10,370 | `add908a387c31886269933366e35bf4a649417d6270cf3acd2fcc1b2d2aaf9c7` | included byte-identically as `docs/signature_algebra.md` |
| `scratch/kobon/n20/obstruction_n20.md` | 16,672 | `e4e223bfc9606624a4c977453a8ddef2cefd5dbfbefd10163d94d23a075a1cfa` | included byte-identically as `proofs/obstruction_n20.md` |
| `scratch/kobon/n20/defect_lemma.md` | 12,371 | `5a1962d2d3aed1af599fe1082b6e2d614cfd4165f516f0ad0f864b9645b26110` | included byte-identically as `proofs/defect_lemma.md` |

## Canonical release contents and archive history

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `math/kobon/release/kobon-2026-08.zip` | live outer container | hash intentionally not pinned inside itself; run `sha256sum` on the archive | canonical release container |
| `math/kobon/release/kobon-2026-08-v1.zip` | 248,593 | `a4dfa51a5c728a49f21fe7c95e8572eaad08b0bca3dd235b13d3d5657c292f63` | referenced only (superseded archive) |
| `math/kobon/release/kobon-2026-08-obsolete-structure.zip` | 248,543 | `052485125b30f3cc31d2ccdd57f9a65d977936ae2aadb5a13f7d61018b3ad3a3` | referenced only (superseded archive) |
| `math/kobon/release/kobon-2026-08-obsolete-structure2.zip` | 248,571 | `bcd90fbcdffd32155d3df701d5568b03e1dd9e2c1ac912ebd0836935d49da2df` | referenced only (superseded archive) |
| `math/kobon/release/kobon-2026-08-obsolete-structure3.zip` | 248,513 | `f29ed2aa4f8a6cbae19a7ebf2c1cd60a60043265b0a132062ce036450d2cb822` | referenced only (superseded archive) |
| `math/kobon/release/kobon-2026-08/README.md` | 3,380 | `d51a6c6c73b7cb76a8968bcab7f5e62b8e9810bbf9dc722ee8936a5d93d129aa` | included as `README.md` |
| `math/kobon/release/kobon-2026-08/discoveries/n11_k32.json` | 3,284 | `c325c10c38a703ed21665afcfbe5159c35510770b4c104850e8e2e58555a23de` | included as `discoveries/n11_k32.json` |
| `math/kobon/release/kobon-2026-08/discoveries/n14_bader53_verified.json` | 1,845 | `2e76304ff6c4c38da88f1b271526a7c00a1b8b50c7168b394dee93880cf27e81` | included as `discoveries/n14_bader53_verified.json` |
| `math/kobon/release/kobon-2026-08/discoveries/n14_best_exact53.json` | 2,689 | `cbbedebdc6eeda51147093b47df34dfe0ffb46f3da5d03fc866c18ebe26fe5f1` | included as `discoveries/n14_best_exact53.json` |
| `math/kobon/release/kobon-2026-08/discoveries/n18_known93.json` | 4,725 | `799a47623b0197109cd70a29b5bb962cff0cd195c34dca2af2a7ea9a9d8cbe2c` | included as `discoveries/n18_known93.json` |
| `math/kobon/release/kobon-2026-08/discoveries/n20_known116.json` | 7,683 | `924e230e9475abffe54271bc418941dd9408ebaa1f7da0677275c30b22d65dc3` | included as `discoveries/n20_known116.json` |
| `math/kobon/release/kobon-2026-08/discoveries/n9_t21_satpipe_rehearsal.json` | 1,342 | `34bd9bc4d3d79af176390813cb56f9d575f5ab9c3fe0060b0910e3f9b055506f` | included as `discoveries/n9_t21_satpipe_rehearsal.json` |
| `math/kobon/release/kobon-2026-08/docs/PROGRESS.md` | 231,594 | `5a5c562ad9761859179e9768691a0bf19c60208cafe5500833e954c50bcf6eef` | included as `docs/PROGRESS.md` |
| `math/kobon/release/kobon-2026-08/docs/batch2_launch.md` | 2,947 | `265635197a56125bfc81ffd5c96b9fcd7c8160f243c0111b553611c3252aa0f6` | included as `docs/batch2_launch.md` |
| `math/kobon/release/kobon-2026-08/docs/litrefs.md` | 16,730 | `12eeacd1f657b8127b88dcf2e1d06e3d43364d13394616879857e4ab6f378239` | included as `docs/litrefs.md` |
| `math/kobon/release/kobon-2026-08/docs/orbit_cube_plan.md` | 12,264 | `69a7b58b9de42a46097ad734251a903265c8b4bfb4006123b9ea897ca173aa1a` | included as `docs/orbit_cube_plan.md` |
| `math/kobon/release/kobon-2026-08/docs/signature_algebra.md` | 10,370 | `add908a387c31886269933366e35bf4a649417d6270cf3acd2fcc1b2d2aaf9c7` | included as `docs/signature_algebra.md` |
| `math/kobon/release/kobon-2026-08/maiorana/LICENSE` | 18,657 | `9ba9550ad48438d0836ddab3da480b3b69ffa0aac7b7878b5a0039e7ab429411` | included as `maiorana/LICENSE` |
| `math/kobon/release/kobon-2026-08/maiorana/SHA256SUMS.txt` | 4,165 | `fc96703b2566e0c6eb6370a911beb436a2e4b01ee26872219de6594c773d0b51` | included as `maiorana/SHA256SUMS.txt` |
| `math/kobon/release/kobon-2026-08/maiorana/sol10_lines_rational.json` | 3,386 | `dc1996c6ee2fe6bbf9d249ef8083cf2ac48fdc564c4ed963d7b3699a8cf2c5de` | included as `maiorana/sol10_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol11_lines_rational.json` | 3,700 | `2310ac2cc0d474d68956e89afec3ed65e164705aa9b8487777736e17d65d8a54` | included as `maiorana/sol11_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol12_lines_rational.json` | 3,374 | `c31c601dd4586302f71075e9f1b302d825d0ecd0500ec6a45d2133f787862e38` | included as `maiorana/sol12_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol13_lines_rational.json` | 3,556 | `6fe8d56686c50309d2e4354791086d7ad7973cc7471a6357f4d8c6757daa4962` | included as `maiorana/sol13_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol14_lines_rational.json` | 3,373 | `d54843bb3bb36d00402bef7cd0ab32bee43ef79de05c7397c6702b8d385398ce` | included as `maiorana/sol14_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol15_lines_rational.json` | 3,424 | `83c4cebf3f2e6baec3748390d55ca4a6e5be5b1017f0b69d3fa44a5aa85d1317` | included as `maiorana/sol15_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol1_lines_rational.json` | 3,541 | `83fc26666d73cd39791018a72aba71caeac4a3f7efb27109ed2955b4540b2523` | included as `maiorana/sol1_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol2_lines_rational.json` | 3,544 | `8d115478b7e1e5324cfb414eb58e261b139c3151f9de941aad6a4df4e0df9524` | included as `maiorana/sol2_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol3_lines_rational.json` | 3,367 | `9212f62b1c72aee7ad874a95fad38575e57b12af331790f141a1a0f2eacc33b8` | included as `maiorana/sol3_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol4_lines_rational.json` | 3,373 | `1bcea050d748df21085599ad6f73b657de828de708e78f92944600e04fd95eaf` | included as `maiorana/sol4_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol5_lines_rational.json` | 3,388 | `2847851c337f13e4d196af1613c766b83ae241ae06de08078f09d92784430142` | included as `maiorana/sol5_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol6_lines_rational.json` | 3,397 | `559f86578d5e77c3b673520174b38670121d1d624cd6e78c2743d0c740be5884` | included as `maiorana/sol6_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol7_lines_rational.json` | 3,369 | `c1129adf88ff74932a93d90541bdad12a5ab99da4209489fc1e933e2dfaa94f9` | included as `maiorana/sol7_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol8_lines_rational.json` | 3,356 | `3c29a2ab722757ceeaa9701f91c4b1a4560043194640f872050d1116fca6620d` | included as `maiorana/sol8_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/maiorana/sol9_lines_rational.json` | 3,538 | `9843d43b6da3d42543f263145b2c5423462cdda0a2afac311f6c95923ef19922` | included as `maiorana/sol9_lines_rational.json` |
| `math/kobon/release/kobon-2026-08/papers/paper_capacity.md` | 60,130 | `f8be81157d75e3f530d194de7e72916134a3d0531a60059b2db652dada3133e7` | included as `papers/paper_capacity.md` |
| `math/kobon/release/kobon-2026-08/papers/paper_kobon_2026-08.md` | 9,973 | `63eeac4be91857da840487c2ef10bcd149c5bbb6d4f5735d8ea75e1c1ebeedef` | included as `papers/paper_kobon_2026-08.md` |
| `math/kobon/release/kobon-2026-08/papers/report.md` | 24,616 | `d6b6e95ac27e1dac478f3dd3018ee7d929cae3862734379eedc0538a14877333` | included as `papers/report.md` |
| `math/kobon/release/kobon-2026-08/proofs/Kgen14_ge54_certificate.md` | 2,970 | `240fd7e18989719923fd7bdc0c328cf96052e86b9c6387110123c0f2cbf7dcfe` | included as `proofs/Kgen14_ge54_certificate.md` |
| `math/kobon/release/kobon-2026-08/proofs/appendix_check_n18.py` | 7,297 | `f423515f9e96a0f2cb495aed54713ae5e05aeb87c25b25450db7705d9754aebc` | included as `proofs/appendix_check_n18.py` |
| `math/kobon/release/kobon-2026-08/proofs/appendix_independent_check.py` | 4,340 | `4292e5ac3b6e17b7daaa6d74d8562695b2c16864a4d429a922c47f574109aa5a` | included as `proofs/appendix_independent_check.py` |
| `math/kobon/release/kobon-2026-08/proofs/defect_lemma.md` | 12,371 | `5a1962d2d3aed1af599fe1082b6e2d614cfd4165f516f0ad0f864b9645b26110` | included as `proofs/defect_lemma.md` |
| `math/kobon/release/kobon-2026-08/proofs/escape_obstruction.md` | 16,862 | `def268d66bcd339ec3a61eaec66dd899106d4fb08c48bd90e2f079386bb7f424` | included as `proofs/escape_obstruction.md` |
| `math/kobon/release/kobon-2026-08/proofs/obstruction_n18.md` | 18,686 | `5da2ddcfaedecab898c520ce6f08da55bde61dd510d75703660944b938ce20fa` | included as `proofs/obstruction_n18.md` |
| `math/kobon/release/kobon-2026-08/proofs/obstruction_n20.md` | 16,672 | `e4e223bfc9606624a4c977453a8ddef2cefd5dbfbefd10163d94d23a075a1cfa` | included as `proofs/obstruction_n20.md` |
| `math/kobon/release/kobon-2026-08/scripts/concurrency_cases.py` | 11,742 | `99d7fca7770cae0dad45294478b2eb14e5b4336815cd332ca6e059e912efc59d` | included as `scripts/concurrency_cases.py` |
| `math/kobon/release/kobon-2026-08/scripts/engine.py` | 43,576 | `2889cace737f26046316395dbf89777d8af51aec65b4291cf2c216a86c7c3029` | included as `scripts/engine.py` |
| `math/kobon/release/kobon-2026-08/verification/maiorana_verify.log` | 520 | `02b0cbe10eac206025de3d68d17711883cde1cc481636d5742e3752eda09e9f6` | included as `verification/maiorana_verify.log` |
| `math/kobon/release/kobon-2026-08/verification/maiorana_verify.py` | 6,220 | `2536ff224256fcb1739ee22291bbf1df39043faee2886d72a86979a781672202` | included as `verification/maiorana_verify.py` |
| `math/kobon/release/kobon-2026-08/verification/sweep_all_verify.log` | 1,601 | `efc7d3e99c8c20812790759c25ef3d917ecec1931e7157646860cec660c23652` | included as `verification/sweep_all_verify.log` |
| `math/kobon/release/kobon-2026-08/verification/sweep_all_verify.py` | 2,960 | `b369be4669be829c5b07ef7668aa2bf6b3ce0d249b8f9468aa26c12f0d0bb401` | included as `verification/sweep_all_verify.py` |

## Core code, copied verifiers, and transcripts

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `math/kobon/engine.py` | 43,576 | `2889cace737f26046316395dbf89777d8af51aec65b4291cf2c216a86c7c3029` | included byte-identically as `scripts/engine.py` |
| `math/kobon/test_engine.py` | 14,022 | `c65bc8c0c0da9098f33f8f5a1467323aa5d2aeddcef57bb209da76cf21651d9b` | referenced only |
| `math/kobon/verification/maiorana_verify.py` | 6,220 | `2536ff224256fcb1739ee22291bbf1df39043faee2886d72a86979a781672202` | included byte-identically as `verification/maiorana_verify.py` |
| `math/kobon/verification/sweep_all_verify.py` | 2,960 | `b369be4669be829c5b07ef7668aa2bf6b3ce0d249b8f9468aa26c12f0d0bb401` | included byte-identically as `verification/sweep_all_verify.py` |
| `math/kobon/verification/max_packing.py` | 5,329 | `c1fba963b52f324e32f46b0f3ba3af05eca3ed7d447c3a9758bbdf2fec9edd59` | referenced only |
| `math/kobon/verification/appendix_check_n18.py` | 7,297 | `f423515f9e96a0f2cb495aed54713ae5e05aeb87c25b25450db7705d9754aebc` | included byte-identically as `proofs/appendix_check_n18.py` |
| `math/kobon/verification/appendix_independent_check.py` | 4,340 | `4292e5ac3b6e17b7daaa6d74d8562695b2c16864a4d429a922c47f574109aa5a` | included byte-identically as `proofs/appendix_independent_check.py` |
| `math/kobon/verification/maiorana_verify.log` | 520 | `02b0cbe10eac206025de3d68d17711883cde1cc481636d5742e3752eda09e9f6` | included byte-identically as `verification/maiorana_verify.log` |
| `math/kobon/verification/sweep_all_verify.log` | 1,601 | `efc7d3e99c8c20812790759c25ef3d917ecec1931e7157646860cec660c23652` | included byte-identically as `verification/sweep_all_verify.log` |
| `scratch/kobon/max_packing.py` | 5,329 | `c1fba963b52f324e32f46b0f3ba3af05eca3ed7d447c3a9758bbdf2fec9edd59` | referenced only |
| `scratch/kobon/n14/maiorana_verify.py` | 6,220 | `2536ff224256fcb1739ee22291bbf1df39043faee2886d72a86979a781672202` | included byte-identically as `verification/maiorana_verify.py` |
| `scratch/kobon/n14/sweep_all_verify.py` | 2,960 | `b369be4669be829c5b07ef7668aa2bf6b3ce0d249b8f9468aa26c12f0d0bb401` | included byte-identically as `verification/sweep_all_verify.py` |
| `scratch/kobon/n14/maiorana_verify.log` | 520 | `02b0cbe10eac206025de3d68d17711883cde1cc481636d5742e3752eda09e9f6` | included byte-identically as `verification/maiorana_verify.log` |
| `scratch/kobon/n14/sweep_all_verify.log` | 1,601 | `efc7d3e99c8c20812790759c25ef3d917ecec1931e7157646860cec660c23652` | included byte-identically as `verification/sweep_all_verify.log` |
| `scratch/kobon/n14/concurrency_cases.py` | 11,742 | `99d7fca7770cae0dad45294478b2eb14e5b4336815cd332ca6e059e912efc59d` | included byte-identically as `scripts/concurrency_cases.py` |
| `scratch/kobon/n14_cube_build.py` | 2,469 | `d1be64f8f73f2e01adde4bbd23e03bca5e4e33bbc00b6ec2f3763839c4274821` | referenced only |
| `scratch/kobon/n18/cube_build.py` | 2,469 | `df92b03b0318a8faf9e06d4971ea4891f06e2a7151cc9d5a57b13878239b7978` | referenced only |
| `scratch/kobon/n20/cube_build.py` | 3,307 | `bfb8926acea9cdd9054721b03089371eefd4493cf6ce7ff57d7e8bd708ce4c4d` | referenced only |
| `scratch/kobon/n18/appendix_check_n18.py` | 7,297 | `f423515f9e96a0f2cb495aed54713ae5e05aeb87c25b25450db7705d9754aebc` | included byte-identically as `proofs/appendix_check_n18.py` |
| `scratch/kobon/n14/escape/appendix_independent_check.py` | 4,340 | `4292e5ac3b6e17b7daaa6d74d8562695b2c16864a4d429a922c47f574109aa5a` | included byte-identically as `proofs/appendix_independent_check.py` |

## Maiorana witnesses and exact packing outputs

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/n14/maiorana/SHA256SUMS.txt` | 4,165 | `fc96703b2566e0c6eb6370a911beb436a2e4b01ee26872219de6594c773d0b51` | included byte-identically as `maiorana/SHA256SUMS.txt` |
| `scratch/kobon/n14/maiorana/sol1_lines_rational.json` | 3,541 | `83fc26666d73cd39791018a72aba71caeac4a3f7efb27109ed2955b4540b2523` | included byte-identically as `maiorana/sol1_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol2_lines_rational.json` | 3,544 | `8d115478b7e1e5324cfb414eb58e261b139c3151f9de941aad6a4df4e0df9524` | included byte-identically as `maiorana/sol2_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol3_lines_rational.json` | 3,367 | `9212f62b1c72aee7ad874a95fad38575e57b12af331790f141a1a0f2eacc33b8` | included byte-identically as `maiorana/sol3_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol4_lines_rational.json` | 3,373 | `1bcea050d748df21085599ad6f73b657de828de708e78f92944600e04fd95eaf` | included byte-identically as `maiorana/sol4_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol5_lines_rational.json` | 3,388 | `2847851c337f13e4d196af1613c766b83ae241ae06de08078f09d92784430142` | included byte-identically as `maiorana/sol5_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol6_lines_rational.json` | 3,397 | `559f86578d5e77c3b673520174b38670121d1d624cd6e78c2743d0c740be5884` | included byte-identically as `maiorana/sol6_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol7_lines_rational.json` | 3,369 | `c1129adf88ff74932a93d90541bdad12a5ab99da4209489fc1e933e2dfaa94f9` | included byte-identically as `maiorana/sol7_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol8_lines_rational.json` | 3,356 | `3c29a2ab722757ceeaa9701f91c4b1a4560043194640f872050d1116fca6620d` | included byte-identically as `maiorana/sol8_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol9_lines_rational.json` | 3,538 | `9843d43b6da3d42543f263145b2c5423462cdda0a2afac311f6c95923ef19922` | included byte-identically as `maiorana/sol9_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol10_lines_rational.json` | 3,386 | `dc1996c6ee2fe6bbf9d249ef8083cf2ac48fdc564c4ed963d7b3699a8cf2c5de` | included byte-identically as `maiorana/sol10_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol11_lines_rational.json` | 3,700 | `2310ac2cc0d474d68956e89afec3ed65e164705aa9b8487777736e17d65d8a54` | included byte-identically as `maiorana/sol11_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol12_lines_rational.json` | 3,374 | `c31c601dd4586302f71075e9f1b302d825d0ecd0500ec6a45d2133f787862e38` | included byte-identically as `maiorana/sol12_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol13_lines_rational.json` | 3,556 | `6fe8d56686c50309d2e4354791086d7ad7973cc7471a6357f4d8c6757daa4962` | included byte-identically as `maiorana/sol13_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol14_lines_rational.json` | 3,373 | `d54843bb3bb36d00402bef7cd0ab32bee43ef79de05c7397c6702b8d385398ce` | included byte-identically as `maiorana/sol14_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol15_lines_rational.json` | 3,424 | `83c4cebf3f2e6baec3748390d55ca4a6e5be5b1017f0b69d3fa44a5aa85d1317` | included byte-identically as `maiorana/sol15_lines_rational.json` |
| `scratch/kobon/n14/maiorana/sol1_lines_rational.maxpacking.json` | 1,620 | `073860dbf6d77c9f324895bf36394a11b749cece01599f18ffc0e8d76a10e4ba` | referenced only |
| `scratch/kobon/n14/maiorana/sol2_lines_rational.maxpacking.json` | 1,621 | `f2f00bfaf1ee6dd448f0b6a2f9ed583491a2a442b1185e3b94311ebcf3806f24` | referenced only |
| `scratch/kobon/n14/maiorana/sol3_lines_rational.maxpacking.json` | 1,620 | `05198f1fb7654a1ef18222be6601cb4a6bb540cf431dd04d81476a49f78fad82` | referenced only |
| `scratch/kobon/n14/maiorana/sol4_lines_rational.maxpacking.json` | 1,617 | `55a5dc9fa0a37ab81c9ce817c04fd9ed2435fae32ac63a87647c9209f3f5cc09` | referenced only |
| `scratch/kobon/n14/maiorana/sol5_lines_rational.maxpacking.json` | 1,617 | `ad93823e35489bed99a3ae3fadf8da255354bc969adb68bc01f2aaf19fe03d96` | referenced only |
| `scratch/kobon/n14/maiorana/sol6_lines_rational.maxpacking.json` | 1,615 | `29b1044f4ebe18530dfb0d529e2d7f8bf9981da78e5fa8acd8b3bba950675465` | referenced only |
| `scratch/kobon/n14/maiorana/sol7_lines_rational.maxpacking.json` | 1,618 | `b819a0d0e25b5da923f2f673010a5b78e03ff91c252defb68398b471194b49c5` | referenced only |
| `scratch/kobon/n14/maiorana/sol8_lines_rational.maxpacking.json` | 1,619 | `b09d4ab2b915c6bc6f168c9b98d8c36ed0d107f4bc7c74fd2eae3fae8c41b843` | referenced only |
| `scratch/kobon/n14/maiorana/sol9_lines_rational.maxpacking.json` | 1,620 | `5b2ec09e2c4dc75eb539bc6bf8a593d3b266f98e53a95c812fe650ad83233456` | referenced only |
| `scratch/kobon/n14/maiorana/sol10_lines_rational.maxpacking.json` | 1,618 | `7ac1c936e6c21a38bd9effed96496dbfe4beb28b45916f02f9b23cd643939d0e` | referenced only |
| `scratch/kobon/n14/maiorana/sol11_lines_rational.maxpacking.json` | 1,619 | `146cc61dcee9cf1b91d655a55807318304783dd333df1a93be816ea4e91e3768` | referenced only |
| `scratch/kobon/n14/maiorana/sol12_lines_rational.maxpacking.json` | 1,616 | `cb8f4503d1378bb4768d64c40f8c131879b400a780d787dded92139799821ba2` | referenced only |
| `scratch/kobon/n14/maiorana/sol13_lines_rational.maxpacking.json` | 1,621 | `2685e32b055af68e88ef0db0d6272fcfff87a0708d69bb519cbc51c0a10dd05f` | referenced only |
| `scratch/kobon/n14/maiorana/sol14_lines_rational.maxpacking.json` | 1,620 | `7ed132b8eaadbbe40523dada61e5f8332ab60bcad4a5e660d042d66ed33a0a9b` | referenced only |
| `scratch/kobon/n14/maiorana/sol15_lines_rational.maxpacking.json` | 1,619 | `150a0f2c3be938b1f6ff52cfb72541b11cbbfb1edb035a4f0a55dce36d40dee6` | referenced only |
| `scratch/kobon/n14/sol1_lines_rational.json` | 3,541 | `83fc26666d73cd39791018a72aba71caeac4a3f7efb27109ed2955b4540b2523` | included byte-identically as `maiorana/sol1_lines_rational.json` |
| `scratch/kobon/n14/sol1_engine_check.json` | 1,927 | `2599147444823c64143dc6f035f6dea447c5f21f67f1505cc09137db9f683736` | referenced only |

## Discovery and fixed-arrangement certificates

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/discoveries/n11_k32.json` | 3,284 | `c325c10c38a703ed21665afcfbe5159c35510770b4c104850e8e2e58555a23de` | included byte-identically as `discoveries/n11_k32.json` |
| `scratch/kobon/discoveries/n9_t21_satpipe_rehearsal.json` | 1,342 | `34bd9bc4d3d79af176390813cb56f9d575f5ab9c3fe0060b0910e3f9b055506f` | included byte-identically as `discoveries/n9_t21_satpipe_rehearsal.json` |
| `scratch/kobon/discoveries/n20_known116.json` | 7,683 | `924e230e9475abffe54271bc418941dd9408ebaa1f7da0677275c30b22d65dc3` | included byte-identically as `discoveries/n20_known116.json` |
| `scratch/kobon/discoveries/n14_incumbents.jsonl` | 5,161 | `d7d3e5d00f928ee6d5363e2bf4084f209dbba01fd455649d9beb48db820d3403` | referenced only |
| `scratch/kobon/discoveries/n14_best_exact53.json` | 2,689 | `cbbedebdc6eeda51147093b47df34dfe0ffb46f3da5d03fc866c18ebe26fe5f1` | included byte-identically as `discoveries/n14_best_exact53.json` |
| `scratch/kobon/discoveries/n14_exact_checked.jsonl` | 3,188 | `35105ac183ef1f1942cb53528bd9c17f747f3386349ecf84ad65c63f004ef1b4` | referenced only |
| `scratch/kobon/discoveries/n14_bader53_verified.json` | 1,845 | `2e76304ff6c4c38da88f1b271526a7c00a1b8b50c7168b394dee93880cf27e81` | included byte-identically as `discoveries/n14_bader53_verified.json` |
| `scratch/kobon/discoveries/n18_known93.json` | 4,725 | `799a47623b0197109cd70a29b5bb962cff0cd195c34dca2af2a7ea9a9d8cbe2c` | included byte-identically as `discoveries/n18_known93.json` |
| `scratch/kobon/discoveries/n20_known116.maxpacking.json` | 3,389 | `81bdf63153804a183ea08a68ba5f0b44ec31ffefe75e589622b9e5bae0f09ef7` | referenced only |
| `scratch/kobon/discoveries/n11_k32.maxpacking.json` | 995 | `6cb44dd990fc96776d3b78dab27dd38e4d8bd112b46816dceefdb9e17dd4cefc` | referenced only |
| `scratch/kobon/discoveries/n18_known93.maxpacking.json` | 2,734 | `38846d512e32c47947938d3f62f23dc47a818c0e13ec8367a9378f4d5a472104` | referenced only |
| `scratch/kobon/discoveries/n14_best_exact53.maxpacking.json` | 1,587 | `336bf3f2e74268ed9e83a65ab0937a01980d77242bec32416c0d7931dc08fbd3` | referenced only |
| `scratch/kobon/discoveries/n14_bader53_verified.maxpacking.json` | 1,591 | `5a10e9fd4601eb7607644218b104eac982b37c5d032933fb5115b6dbfe0015c2` | referenced only |

## Current n=14, n=18, and n=20 CNFs, logs, and live proof streams

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/n14_monolith_t54.cnf` | 240,994,798 | hash omitted (large) | referenced only |
| `scratch/kobon/n14_cube_q0_t54.cnf` | 198,349,702 | hash omitted (large) | referenced only |
| `scratch/kobon/n14_cube_q1_t54.cnf` | 197,994,509 | hash omitted (large) | referenced only |
| `scratch/kobon/n14_cube_q2_t54.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14_cube_q3paired_t54.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14_cube_q3triad_t54.cnf` | 197,282,896 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/n14_monolith_t55.cnf` | 240,627,409 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/n14_monolith_t56.cnf` | 240,259,861 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/n14c-q0tp1.cnf` | 240,995,268 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/n14c-q0tp2s.cnf` | 240,995,273 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/n14c-q1tp1d.cnf` | 240,995,271 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/n14c-q1tp1s.cnf` | 240,995,270 | hash omitted (large) | referenced only |
| `scratch/kobon/n18/n18_cube_q0_t94.cnf` | 1,063,540,279 | hash omitted (large) | referenced only |
| `scratch/kobon/n18/n18_cube_q1_t94.cnf` | 1,062,343,488 | hash omitted (large) | referenced only |
| `scratch/kobon/n18/n18_cube_q2_t94.cnf` | 1,061,146,329 | hash omitted (large) | referenced only |
| `scratch/kobon/n18/n18_cube_q3paired_t94.cnf` | 1,059,948,802 | hash omitted (large) | referenced only |
| `scratch/kobon/n18/n18_cube_q3triad_t94.cnf` | 1,059,948,532 | hash omitted (large) | referenced only |
| `scratch/kobon/n20/n20_cube_q0_t117.cnf` | 2,131,260,256 | hash omitted (large) | referenced only |
| `scratch/kobon/n20/n20_cube_q1_t117.cnf` | 2,129,377,829 | hash omitted (large) | referenced only |
| `scratch/kobon/n20/n20_cube_q2_t117.cnf` | 2,127,495,034 | hash omitted (large) | referenced only |
| `scratch/kobon/n20/n20_cube_q3paired_t117.cnf` | 2,125,611,871 | hash omitted (large) | referenced only |
| `scratch/kobon/n20/n20_cube_q3triad_t117.cnf` | 2,125,611,601 | hash omitted (large) | referenced only |
| `scratch/kobon/n20/n20_cube_q4four_t117.cnf` | 2,123,728,340 | hash omitted (large) | referenced only |
| `scratch/kobon/n20/n20_cube_q4pairtriad_t117.cnf` | 2,123,728,070 | hash omitted (large) | referenced only |
| `scratch/kobon/n14_monolith_t54.kissat0.log` | 2,235,847 | `bdacdc0f1b3ab804f770a12074091efcfdcd7023bde9151f3ace19919d4dc0f4` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_monolith_t54.kissat1.log` | 187,952 | `47597e924dd6c890287d4dd756828c8d5218ddaa98635f485bb9118201702024` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_monolith_t55.kissat0.log` | 1,295,060 | `3a9de9d55fffe367b95b9998b4fd6b567c774cec55ddd77f9a16d213b0b1becb` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_monolith_t56.kissat0.log` | 102,669 | `22e69adcae83301b65ea87a70a36d442381f751d375bcf6b217fced76933c4dd` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_cube_q0_t54.kissat0.log` | 2,213,538 | `603714be90689227b8e8100be03e21da06f92fce231d3f875152b12b6a93a9bf` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_cube_q1_t54.kissat0.log` | 2,305,893 | `96a58dc7bcdb6be6e63ae0f9b27dcd2adfd5f1572866c28106f6228bf244680d` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_cube_q2_t54.kissat0.log` | 2,227,819 | `9d3f7aeda0a2e0a9c52cb62a4890d6de85ad16c0d1eeb4958df881edbc789348` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_cube_q3paired_t54.kissat0.log` | 280,242 | `72bf0ac9147742e84f81fddba83b48f9062ea25c4a070765b04f98ba81bd8f12` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_cube_q3triad_t54.kissat0.log` | 2,649,159 | `56d439c1a2db9bafcfcf00e386baa3c8dc3ec14259d7daeb2b7769b3d2662510` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/n14c-q0tp1.cnf.kissat.log` | 1,924,810 | `5a059b600ed8f9f35a0b3877dd7dfbc8c26c78f9418eec95514fff6827d878ff` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/n14c-q0tp2s.kissat.log` | 1,258,107 | `760f8a05cee38f3ba86e8da945f2ce4ae1c054ed3e586b3d5c423bc54afad760` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/n14c-q1tp1d.cnf.kissat.log` | 1,849,933 | `ed2a6401e55b6d51ca1974377b2a73a81d7d3402e8df795169b7cbdb44cee633` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/n14c-q1tp1s.cnf.kissat.log` | 80,030 | `c9b9b02da6906209624cc16c3bfa796db936a09a743abd8966d5d2b9a5e00d4e` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n18/n18_cube_q0_t94.kissat0.log` | 813,923 | `0cc1b3d1c3a56c5a65a263fd153b0ad55e21902c44528c9e667dbdd65ab0dc50` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n18/n18_cube_q2_t94.kissat0.log` | 889,335 | `5fb780596060cccf1a53958222f119f186413ec0f60e37021c453889a8f9f941` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n18/q1_cadical.log` | 129 | `183624e69ca1e8710faacf3a294439ea70a52ef5c50ee4f18863b146cc19f4eb` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/n14c-q0tp1.cnf.drat` | 75,473,354,752 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/n14c-q1tp1d.cnf.drat` | 81,948,311,552 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/n14c-q1tp1s.cnf.drat` | 464,519,168 | hash omitted (large) | referenced only (live/incomplete at snapshot) |

## n=14 generated position cubes and launched logs

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/n14/q2pos/q2pos-000.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-001.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-002.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-003.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-004.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-005.cnf` | 197,638,989 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-006.cnf` | 197,638,991 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-007.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-008.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-009.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-010.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-011.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-012.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-013.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-014.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-015.cnf` | 197,638,989 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-016.cnf` | 197,638,991 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-017.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-018.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-019.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-020.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-021.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-022.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-023.cnf` | 197,638,989 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-024.cnf` | 197,638,991 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-025.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-026.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-027.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-028.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-029.cnf` | 197,638,989 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-030.cnf` | 197,638,991 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-031.cnf` | 197,638,992 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-032.cnf` | 197,638,988 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-033.cnf` | 197,638,989 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-034.cnf` | 197,638,991 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-035.cnf` | 197,638,989 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-000.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-001.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-003.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-004.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-005.cnf` | 197,283,140 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-006.cnf` | 197,283,142 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-007.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-008.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-009.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-010.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-011.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-012.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-013.cnf` | 197,283,140 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-014.cnf` | 197,283,142 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-015.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-016.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-017.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-018.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-019.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-020.cnf` | 197,283,140 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-021.cnf` | 197,283,142 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-022.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-023.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-024.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-025.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-026.cnf` | 197,283,140 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-027.cnf` | 197,283,142 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-028.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-029.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-030.cnf` | 197,283,139 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-031.cnf` | 197,283,140 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-032.cnf` | 197,283,142 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q3pos/q3pos-033.cnf` | 197,283,143 | hash omitted (large) | referenced only |
| `scratch/kobon/n14/q2pos/q2pos-001.kissat.log` | 1,754,192 | `3ae2d8b6ba060378118a9c7e6c0debf5f287cb16f3a547b30fbb8b1a02c31e57` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/q2pos/q2pos-002.kissat.log` | 1,539,558 | `069a9500f35cab52955d5c586278d919d34fa2413a60d4f0ee3156f6c7bfd5b1` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/q3pos/q3pos-000.kissat.log` | 355,039 | `219b95d86ff7178b45d703ec209b7ae6f8b09a65ee6229797da1cf9753acbaa4` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/q3pos/q3pos-001.kissat.log` | 275,027 | `5283677c9eaafd52838c4901789a8fd365dad89c8193db6e2830210fed8b8781` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/q2pos_index.json` | 15,291 | `25503ffee66d40a54d1942550eb7dafbc4cacd259d8b4d3c0615fcc98e257837` | referenced only |
| `scratch/kobon/n14/q3pos_index.json` | 34,127 | `c5509239be6e86c10539225413889ac98bab6be0cb9ba7d97c67d47158b5c1a8` | referenced only |
| `scratch/kobon/n14/q2_positions.py` | 8,792 | `7b66ce4080904a1d37b1cccef328a232f1c575028ef2ddcea171972b07622248` | referenced only |
| `scratch/kobon/n14/q3_positions.py` | 12,027 | `2c6b7a5da6add0e05c0f2348196949fa689d2d9d18791424cf7be77c7974dff3` | referenced only |
| `scratch/kobon/n14/q3_radius_probe.py` | 4,706 | `9396db2a6fe3d26d0678fa01731a1980ddf4b48ef79776c3d1c73f821ff2042d` | referenced only |

## n=11 witness, analytic microchecks, and live monolith

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/n11_monolith_t33.cnf` | 44,481,369 | `dac2849cf3ed19fa8e5e17224b9a9d594206bfb91d12c1b8f32017411be77491` | referenced only |
| `scratch/kobon/n11_monolith_t33.drat` | 94,648,664,064 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n11_monolith_t33.kissat.log` | 1,631,571 | `ecf68f19441624a4156b14c0649807be70575f2e4007755ceec5bcea53cdda23` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n11/flower_n11_verify.py` | 3,100 | `e645eb42d396ad3e72850b833672e9cc47dae5084869e95a01f1d08aa60ff514` | referenced only |
| `scratch/kobon/n11/flower_n11_config.json` | 380 | `6bb4c668cb784d2050255020e6863bd1cb93826e6a0c9e0d6ab59195c062b748` | referenced only |
| `scratch/kobon/n11/flower_experiment.py` | 2,305 | `894f17ca3b8e2500ec21f0adf2bd0ccd4942b5164921eb826ec53de2f378d224` | referenced only |
| `scratch/kobon/n11/face_tracer.py` | 4,871 | `5f82ade062ef880ed843f7671e2780d7e74172947ceb03fa520cca2a937ec2d3` | referenced only |
| `scratch/kobon/n11/cube_q0.py` | 2,128 | `543b3196765ca44d9a11a1ea211e6f406e9bc1d36503d14ac33733c4aef18517` | referenced only |

## n=14 obstruction support artifacts

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/n14/escape/multigraph_feasibility.json` | 35,698 | `9b8e70a2e12c6ce188df9aecf8e797cb01384b282e11653f1613898e4204444d` | referenced only |
| `scratch/kobon/n14/escape/escape_theorem_verify.py` | 5,455 | `58686f175be2ebda435ee0f6217e9708f0d3856ea6eb693c77fd0d76c917ec0f` | referenced only |
| `scratch/kobon/n14/escape/bader_snapped_report.json` | 1,161 | `4a0e20d1cb5d6aa5baa690ebb47f5abdfa90cc21618ea0d794bd24af8b017617` | referenced only |
| `scratch/kobon/n14/escape/bader_stored_report.json` | 1,660 | `577b22c2c78f6a29b6083b22dccbea91f2648cebb9ce92011912090b53cc8eec` | referenced only |
| `scratch/kobon/n14/escape/random_stress.json` | 348 | `6a3c491f3b75aab702692afed0e79dd756d86b9eefd38d1ddd920a84eaa8c184` | referenced only |
| `scratch/kobon/n14/escape/escape_random_search.py` | 2,123 | `b1430744355a636372b670cea3ee5ce7f02b631dba704624a256f2e204eb4346` | referenced only |
| `scratch/kobon/n14/escape/escape_multigraph.py` | 2,466 | `3be6636e7d60e6bd11faa015a542585ad37d630345aefb7c174254ec02b34dbf` | referenced only |
| `scratch/kobon/n14/escape/n14_q3_pattern_search.json` | 36,512 | `e16685b56462db5ee380b3f09aab26e756e0efd9d7ce58cfaae0496a657f8dfe` | referenced only |
| `scratch/kobon/n14/escape/escape_n14_q3_search.py` | 2,563 | `913c895eb814074d01fa69aeee3bb0592a512eefbaf3e385cc9f4eff46d8c089` | referenced only |
| `scratch/kobon/n14/escape/small_cases2.json` | 2,453 | `4dffc2a3bd65dc68e7d230bdcd5ce90b777e0bfa462d9f96c726e9f2e62ee343` | referenced only |
| `scratch/kobon/n14/escape/escape_small_cases.py` | 3,242 | `b9310ba550fe8fc7efcac923872fede68e220251002f2e59a6fe3b0184b4215c` | referenced only |
| `scratch/kobon/n14/escape/small_cases.json` | 4,282 | `7f9958df9750e92e7ba0fcdf389bcd69c6acbf39168b35cc94b3f185d19225fe` | referenced only |
| `scratch/kobon/n14/escape/escape_lib.py` | 7,652 | `af6bda107ca498521b95fc1eeccbe944351e1f78aec7b1d25a94120239f638e7` | referenced only |
| `scratch/kobon/n14/escape/escape_bader_check.py` | 4,838 | `c9731ee1f7ba4889d023212c9928c5e59d44b11cfe69be125d61d3f4ef6695be` | referenced only |

## Savchuk probes and n=9 rehearsal

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/n14_savchuk_global6.cnf` | 6,866,365 | `3df185fc94b7ed90eaf272bd167c3880ad831796261c05d54385853d044659bb` | referenced only |
| `scratch/kobon/n14_savchuk_global6.raw.cnf` | 35,766,302 | `8ed09738d437805e974e686b698b5aabe8b4a9f8ee20271f6ee83ef8e9fd87fc` | referenced only |
| `scratch/kobon/n14_savchuk_global6_exact_sym.cnf` | 9,766,870 | `2fcbea3bd5190d9e71686bb455b33122436923bbfeb6cf07e0755d21799b38ad` | referenced only |
| `scratch/kobon/n14_savchuk_global6.kissat.log` | 184,147 | `01b38dcd6ba187b8f66cf1d2d71fa497d58223b7d7bb8c9a19df3d7d2102550e` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_savchuk_global6_exact_sym.kissat.log` | 171,053 | `cba35de4715022b0f995bbf6e842f2a6abf466eb56cc7d8c7697903b0ed8efda` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14_savchuk_global6-s1.kissat.log` | 96,833 | `59e83498ab5509ea817b019472a1778f76ab211f3d44e47a0cc82489e3cfa16d` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n14/work/probe_results.log` | 5,899 | `ac85f5d95c52b200c971af55c57c30e09b8d6bb31f3d820aeeb7686246f1ae0b` | referenced only |
| `scratch/kobon/n14/work/stack_validation.json` | 680 | `1c4b5d17b5595b963d6fbf07be3e890565ccfffb60a25f92c510d9d6873ac56b` | referenced only |
| `scratch/kobon/n14/work/validate_stack.py` | 5,853 | `17a2fc4a498aa816f502bccdb3cb971b729330a42cc23fc7c98a351869210d9f` | referenced only |
| `scratch/kobon/n14/work/straighten_probe.py` | 4,074 | `553f3bea5b340bb9e54b0eb79bddb5708206e3fb442a39ce877f53e53fb5c5e8` | referenced only |
| `scratch/kobon/n9_rehearsal/n9_t21_recipe.cnf` | 10,055,181 | `aa67867462cd20ba9d2ef2a70a1d2c8f89943c54c95248599661e4bdf30b444b` | referenced only |
| `scratch/kobon/n9_rehearsal/n9_t21_recipe.cnf.model` | 405,552 | `c0954844a5f78767eb107a9730d6306b062f4207e501ac8b51e75936d5537bbf` | referenced only |
| `scratch/kobon/n9_rehearsal/n9_t21.kissat.stdout` | 567,531 | `43f1ddaef8bae2ac0f89246ca87e91eeb3234d3b37989c8b7fba72ac03946253` | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/n9_rehearsal/n9_t21_chB.json` | 1,336 | `41d5363acf71a74f27b4c9c2d82f2a57efc73be5b26089179735c4c6a58af691` | referenced only |
| `scratch/kobon/n9_rehearsal/n9_t21_chC.json` | 1,345 | `de8a17350d74e328187bf79a636dbef3bd3305f1022b962c6777a0d848a1ab52` | referenced only |
| `scratch/kobon/n9_rehearsal/rehearsal_report.json` | 2,758 | `e09d129e5b4eaac0487cfb8d273d3dfaf5616601a7c155f2b2648928f38bb93b` | referenced only |

## n=18/n=20 reconstruction and diagnostic artifacts

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/n20/max_disjunct_optimum.py` | 3,329 | `eee1643e33ae604385bf1e1fe133c51f3cf342097bf35cd3121e277c1335ffd6` | referenced only |
| `scratch/kobon/n20/known116_reconstruct.py` | 9,817 | `05ba60b3c7520dfd1686ad1aeb4054eecf2bbfd2bb7bce77eef642963b3a6e0a` | referenced only |
| `scratch/kobon/n20/known116_reconstruct.log` | 1,004 | `0fae5d2a386155a7b5a35afc07b657b3c58529f6e7e0a54b443425e0f2732348` | referenced only |
| `scratch/kobon/n20/explore1.log` | 847 | `b835e97dcd78d7686598943be75d853aa5b008031d2326a5b56a34c46930e9e0` | referenced only |
| `scratch/kobon/n20/explore2.log` | 207 | `be366bb7df85fc5b7fa8f8d3e2ca0e0ecfd2d5ad4a5c1e86d03b96993a89f985` | referenced only |
| `scratch/kobon/n20/explore3.log` | 765 | `7bbcf8a1523c29c1744527b55a9f189e785de317e127ac4fbe86f65f38b7d8d9` | referenced only |
| `scratch/kobon/n20/explore4.log` | 2,954 | `8871740427d525780bb1fdfe33f87e4f05806fff484328e9a2cc91e4cf94d352` | referenced only |
| `scratch/kobon/n20/solve_cadical.py` | 1,294 | `0717a364ba5361fae2b98b8895ae7f5079ca2b4bfd763efa8d55f438089439e3` | referenced only |
| `scratch/kobon/n18/run_cadical.py` | 984 | `e6dd14c96a58162ff1c78a90ccf43fcde9a4baca39e9d5d6d0658a8f471d2fb1` | referenced only |
| `scratch/kobon/n18/certify_known93.py` | 3,668 | `8c39b82de3ef0c309d24cdb883cbff7a9a31cef025134bde525e82130c744908` | referenced only |

## n=12 retained lower certificate and verification transcripts

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/n12/n12_lower_certificate.json` | 1,590 | `41f2d3460ad72b847e00ef15266e89e2b3c969c9a94ee7f85126b10cf8e97293` | referenced only |
| `scratch/kobon/n12/kobon_n12_t39_proof.cnf` | 64,422,412 | hash omitted (large) | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c0.log` | 237 | `2d93ef7893c1ab4ddc4ea963a7bb061c402aa8c2c63f6ac0a39a44a8c5a2e5ba` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c1.log` | 237 | `5eb69ec60cfb14a985e32427110c8b981b6e595fc21903d2ed0ffea4f7d003c5` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c2.log` | 237 | `1577803cf3ad55f7343f1964b646bb404bc2dd80dbd3f2977980f60b5043d9c3` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c3.log` | 237 | `61001f76558ee54fa997786f3b51148fc50e6b6556b206b17da5836e22482e8d` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c4.log` | 237 | `3269413edc7c71f61ed2a2832c8276cb04684ba52747fa26ae54459d13b592bc` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c5.log` | 237 | `d7573ed98488a63617e06c6bc031eff0a4898a607783cef031313ed25f1e6767` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c6.log` | 237 | `06d3f34b30e750f98fec798e32795c7d1d1cdfc1df8b224a4914c8bb0858a7c0` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c7.log` | 237 | `ad4b669d273e593ecee389906279e29e4a734c45e9484802e1325d51efebdb9f` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c8.log` | 237 | `9e2dfcf286adc04967a98c1264fc1b9116d382eb978f00728d99a0ff9035eaaa` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c9.log` | 237 | `a91c95e888b5763ceb5ee97a654144e1b3ca69fcc8ae1ecdda91e410d0692003` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c10.log` | 237 | `c159624e62879588d49b65072f095a83086d461683d844136f4a56f58a1f51ef` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/c11.log` | 237 | `3ce465193a319a02cabc740e2a45da99c88ea1f41c89e458edde214f7dde8e81` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/simple.log` | 237 | `70ef2bd150d838137c6fe7bcae4c5db9761327b1a0c127b679bfce2086a45966` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/par_noc.log` | 237 | `713e6054662919ab63fd4ca398caf0b5a94202d6872f7881520a7a0ce279dd18` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/split_par_conc_0.log` | 237 | `55248b276cecad15db598f5171347539b8d3dad8045e8b3b36fc88cef8cd3a6c` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/split_par_conc_1.log` | 237 | `775a3ac97136556c87498eaab7f3ff8d93da47f2246f29330d0c17c5577732b0` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/split_par_conc_2.log` | 237 | `fbef648a819c9efe0e4e4cbf8b9430bd600f2aa3a806a3bca9e32fdacea2b445` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/split_par_conc_3.log` | 237 | `1595bb9eaf6a282a83230596b7af9bd83da5de60aa67e656f81cbdafaa84bf6c` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/split_par_conc_4.log` | 237 | `92a0a79e8d7fa9c0e7f14618e04be195a744d55a73c870ad05847702de5c9e02` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/split_par_conc_5.log` | 237 | `c46306138504df46b19142ebdd4bf0b42203322e8942ae0630de3764005a6900` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/split_par_conc_6.log` | 237 | `7b62ab974edb9f4ea6932f6d11d8cb6d09936861e7f4851f913f180fffe3097a` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/split_par_conc_7.log` | 237 | `7a585cfbf0009ad04cbd18082358e9be891163ac08771e0b85adb4161c0112cb` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_0.log` | 235 | `e026dcd144c7bd6ddd0404dce444cfdc838ece85b7ccb7b6cfc2c547e4734c34` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_0_s2.log` | 235 | `b116414b168f4843c66088f573ade545ed6c90f284ca12a60140f53544bafb26` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_1.log` | 431 | `c7747dd39221ae8da1b9d99e0d6c1314bab954ada9a63007bde0b127d40f3e68` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_2.log` | 235 | `8b11b153d9b6be5ce9d77580d82f11fb283944589d48be8ca463503f4f08c965` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_2_s2.log` | 235 | `8c9484b30baaa9e3190cdf2d620cb1b2e533013987f80399bcebbebea816e1c8` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_3.log` | 431 | `c4e7973f046b59338510374669aa50d4c8c91715b7f80aaea1fc212d265820ef` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_4.log` | 235 | `9812c4e9e9b26b7c2751949379a6a2b301c5f8943ad527a1f8e47698a35169fd` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_5.log` | 431 | `872641ab5c51857779708c5aee0ae174ef56eb7388aef7aa08a34687190b05f6` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_6.log` | 235 | `8573f4a16faa97d82dd9b3ba4468580739a4d8ee2868865ed287800ca55e0067` | referenced only |
| `scratch/kobon/n12/verify_attempt1_20260816/deep_par_conc_0_7.log` | 431 | `3948d66ef224cea89886965d2555922fc6db8591df9b80d3b7af3e08ba9ace5c` | referenced only |
| `scratch/kobon/n12/verify_xz/smoke_xz.log` | 230 | `20fc16a625eb09b9dd5852dcacb149e3a7e7d7163bc5b38afd02794284a4a659` | referenced only |
| `scratch/kobon/n12/verify_xz/seed2_smoke.log` | 231 | `dd7c6de2138a52586d2e8528b85ffa469a99d8daee762360ca345d4496550213` | referenced only |
| `scratch/kobon/n12/verify_xz/parnoc_1_a2.log` | 429 | `5afa85290fb2ff267465b82213601ccde5745c0e4fde3ff323f57c5e57079ec6` | referenced only |
| `scratch/kobon/n12/verify_xz/parnoc_3_a2.log` | 429 | `c0834ca3b4eef23ae8063d94c6090247f99b7ec54775d89991f93e4a84e55088` | referenced only |
| `scratch/kobon/n12/verify_xz/parnoc_5_a2.log` | 430 | `72009cb5981635bc4bbdd5c7af4b9df2a56a9128a9834fdc7f135420de13f3ad` | referenced only |
| `scratch/kobon/n12/verify_xz/parnoc_7_a2.log` | 430 | `e9e95bac5a22a5ba8f73d81a1414e8d2dfee6e8489428aded2d4c75b8bb35431` | referenced only |

## Certified n=10 parent/cube and replacement split evidence

| Path | Size (bytes) | SHA-256 | Canonical release zip |
|---|---:|---|---|
| `scratch/kobon/kobon_n10_t26_proof.cnf` | 17,656,901 | `1605080734b32edb31d6a63145de23dc47348544a25617fb206d37750e961a7f` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_simple.cnf` | 17,738,152 | `e29e4cc8b2a03ec1b90694f68b4954ec96071590a1b40991d7936d4259af94a0` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_par_noc.cnf` | 17,736,727 | `7b76ccaac640960f7c04743cbb93370e2fb3729ad028866b86191363a0b4f054` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_par_conc.cnf` | 17,657,332 | `2be61b00821d9dbaea17d68c8869eda024ba8b5f1264510d3bf498ee3ec53203` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_c0.cnf` | 17,657,123 | `031bcbfd4da94536e460795f0cfef9c89ee02fc4a9830435f05a89bfd6c8e48b` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_c1.cnf` | 17,657,123 | `67615ebcfb4036fa871383f7760707921f75ce721701e94324a203517efe98a1` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_c2.cnf` | 17,657,124 | `5909ada386099d0f25bf9cb530caa3f062268327a808da3b32434397bb8f1e8f` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_c3.cnf` | 17,657,124 | `c0ed5b53c66bdea5d07880d53c4228b0fee8351b9dfb98fe586a0dfc3328be9e` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_c4.cnf` | 17,657,124 | `bd6c3e5ddf36e38e69f49ac60c2e65aa1f9e58413b4fd2b4676f0f2c88ea603b` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_c5.cnf` | 17,657,124 | `613761915887c91299f863c912370fd29683e68fc3861b5fe78207c4382b4559` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_c6.cnf` | 17,657,124 | `f2af79254904be68dc6d98c31e9bd9639d116a0bdf1e39d81776d61f1736f63b` | referenced only |
| `scratch/kobon/kobon_n10_t26_proof_c7.cnf` | 17,657,124 | `99e78c3ba2a019dae979df0c9c4df73bb737d27ec0b6d1649b15dcbf7a6fddda` | referenced only |
| `scratch/kobon/proof_n10_t26_simple.drat` | 3,783,608,681 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_par_noc.drat` | 1,777,267,327 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_par_conc.drat` | 10,075,766,784 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_c0.drat` | 8,831,107,072 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_c1.drat` | 5,462,532,361 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_c2.drat` | 6,700,287,805 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_c3.drat` | 5,818,190,931 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_c4.drat` | 5,908,208,736 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_c5.drat` | 7,918,616,844 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_c6.drat` | 7,883,047,474 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon/proof_n10_t26_c7.drat` | 7,298,778,122 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/results.tsv` | 292 | `373687af66f44d0a8b701c1cdc99e31d25ddcdffe466a6e3d6dd0efd4193b115` | referenced only |
| `scratch/kobon-audit/results_verified.tsv` | 308 | `2759863631206a1eb23b87122f777d22ff6b26b61bc139656efea470c22ac690` | referenced only |
| `scratch/kobon-audit/splits.tsv` | 539 | `496ed9b1521ebb8c8ef6702ee3703083eab7587086cd6ec20a9219d7ce0fa3fb` | referenced only |
| `scratch/kobon-audit/logs/simple.log` | 1,260 | `3d3905daf628eeea909f4ab2b4ad2b9d7d61f4d35e87553c9bb33905989667bd` | referenced only |
| `scratch/kobon-audit/logs/par_noc.log` | 1,259 | `28157dc987c638c45579cbd6f98f08cd6dad0800c401bb4e965202972d164d9e` | referenced only |
| `scratch/kobon-audit/logs/c1.log` | 1,262 | `4d66e620e881e3dc44f0a76354b12d9d82de094a1067befd2e0b63376fd04298` | referenced only |
| `scratch/kobon-audit/logs/c2.log` | 1,263 | `bb939cf4d485f25a407fc18ee98340dbe7c200521b206e412987482fab60912a` | referenced only |
| `scratch/kobon-audit/logs/c3.log` | 1,262 | `d97940da1cfda483a6d7a7d08667baff70335c134f0049d58f0f2cc9c8485b35` | referenced only |
| `scratch/kobon-audit/logs/c4.log` | 1,262 | `b417019de33421e5ddf2d64f67a7268b45066dacec210e37b277513e9f3e76d0` | referenced only |
| `scratch/kobon-audit/logs/c5.log` | 1,263 | `87d27ffa2cff668305b038a2dc829f4a230e3560b6e269d0395f611f30bc8c0d` | referenced only |
| `scratch/kobon-audit/logs/c6.log` | 1,263 | `08f97e76fc4762b9ec887789d1f9efda12a69a6af37e84f4010c5ae9e7b45dcf` | referenced only |
| `scratch/kobon-audit/logs/c7.log` | 1,263 | `6f5150ceeed795c31aab4f48bba65d789701511b9f37df84a67fbaddca9a69e3` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_c0_0.cnf` | 17,657,144 | `6055dab01c0d55de572af05aae0e3984221c684fa152873a76a598f7211ee4ec` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_c0_1.cnf` | 17,657,143 | `98f4b1278a07799aa4aed035f423828a3ce5fa975bced2e15a3a410f7d0d4a99` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_c0_2.cnf` | 17,657,143 | `388f6777aa8bbcfd98154cf6a6db2ca642c2b06cc208e40843afc0ee4413cda7` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_c0_3.cnf` | 17,657,142 | `5f271fb96fa71b298f5d24d7ff92cb7b7aba452adbddb522526d2f203969b403` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_c0_4.cnf` | 17,657,143 | `03712ddcaf36b9fe692c232722d3c9eb637c78ce5ffca896d3c7789c55e54d7b` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_c0_5.cnf` | 17,657,142 | `2b80ee0675dea26d53a0842e2d424ca4e5528cfd32e85f4304a6da093f6c2a96` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_c0_6.cnf` | 17,657,142 | `59f4b214d04e93dba9cff88e920713135e6eead9fd399077128c53b076be1925` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_c0_7.cnf` | 17,657,141 | `89bb790f1c75531d73d4bcd1abfac973dd587248382c0f6d07887d7ef95e16a7` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_par_conc_0.cnf` | 17,657,353 | `9e2ab85d9cd8bc316e081c1f0abfd693a1474146bfba6765dc3441746eb619c0` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_par_conc_1.cnf` | 17,657,352 | `71f22aef0b44e3dcd71af45f053362cfd146f910a1904b20ee04eb6522b7c5d9` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_par_conc_2.cnf` | 17,657,352 | `dbf401e3a40bb7e18bf4b867a65d7583c2d4554f542f8e4c189a2dfe80189e97` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_par_conc_3.cnf` | 17,657,351 | `8db6fdc24c77a1d8c0fad63adf82e2a107568075c99f880a2fb4c3760344d4ec` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_par_conc_4.cnf` | 17,657,352 | `b5eba915cd6c3bba524cd0f23029a0543e04058d7e6e43ce6b29c885e1ffe1b6` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_par_conc_5.cnf` | 17,657,351 | `9c6571be5b133c9f30d933a5b4a4ed90a64d5ef2a761523b22d1524b3fbe92df` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_par_conc_6.cnf` | 17,657,351 | `08820ae996ca8d599e480714b0feaf0b5edd20d5dfad7791361997e8a71d9ef8` | referenced only |
| `scratch/kobon-audit/split/kobon_n10_t26_split_par_conc_7.cnf` | 17,657,350 | `c6e63709f2997c0de8caba6760a3203a6c3b773694facafbc796b53aedba9191` | referenced only |
| `scratch/kobon-audit/resolve/split_c0_0.drat` | 1,509,258,056 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_c0_1.drat` | 1,494,763,651 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_c0_2.drat` | 1,096,736,107 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_c0_3.drat` | 1,316,453,716 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_c0_4.drat` | 1,432,999,823 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_c0_5.drat` | 960,574,339 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_c0_6.drat` | 1,422,207,598 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_c0_7.drat` | 1,555,035,128 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_par_conc_0.drat` | 6,802,562,475 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_par_conc_1.drat` | 3,871,676,844 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_par_conc_2.drat` | 5,978,346,362 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_par_conc_3.drat` | 5,425,866,492 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_par_conc_4.drat` | 4,464,633,595 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_par_conc_5.drat` | 2,090,456,056 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_par_conc_6.drat` | 4,243,039,771 | hash omitted (large) | referenced only (live/incomplete at snapshot) |
| `scratch/kobon-audit/resolve/split_par_conc_7.drat` | 3,528,039,283 | hash omitted (large) | referenced only (live/incomplete at snapshot) |

## Requested paths not present

The generator checked the following expected names but did not inventory them because no file existed at the snapshot:

- `scratch/kobon/n14/q3pos/q3pos-002.cnf`
