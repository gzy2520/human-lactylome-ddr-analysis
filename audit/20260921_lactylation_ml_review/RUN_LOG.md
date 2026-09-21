# RUN_LOG — 20260921 lactylation ML re-review (UTC 2026-09-21)
Worktree: /Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912 (branch research/kla31-rnaseq-20260912)
Seed 25; R 4.4.3; data.table/digest/ggplot2 (+ROCR present but unused; base glm used).

- Preflight PASS (163/17/174; 28 refs/1898 samples/31 mappings). fork fd3ace5 used as audit reference only.
- 01a G1G2: G1 holds (no local intensity values); G2 recompute conn=19 (matches fd3ace5).
- 01b G3G4: folds_BRef_leaveConnOut.csv (19 folds, all test N>=20 & +/-≥5 → 19/19 ok, min test pos=5); G4 regroup-manageable.
- 01c G5G6G7: G5 CORRECTED (28 matrices local + 128-gene table); G6 softened (ρ=0.189 p=0.308); G7 alternative implemented.
- 02 BRef task: BRef_task_table.csv.gz (Ref-restricted DDR; rows=6744 pos=2708 rate=0.402); feature_sets_frozen.csv (364/128/48; note 364 vs 371: 7 panel Ensembl absent from per-ref matrices — used complete-case per row, NA=0 rows).
- 03a features: BRef_features.csv.gz (533 genes x 28 refs medians; NA target rows=0).
- 03b train: BRef_oos_preds.csv.gz + BRef_pooled_perf.csv + BRef_perfold_perf.csv. NOTE pooled-M0-ROC artifact (0.32) from cross-fold constant baselines → report fold-mean (M0=0.500) instead; pooled file retained with artifact documented here.
- 04 figures+sensitivity: fig_BRef_perfold_ROC.{png,pdf}+sourcedata; fig_BRef_calibration_M1.{png,pdf}+sourcedata; sensitivity_dropOneConn_M1vsM0.csv (Δ=0.139, range 0.129–0.151).
- Results: fold-mean ROC M0 0.500 / M1 0.639 / M2 0.637 / M3 0.634; PR M0 0.345 / M1 0.492; Brier M0 0.2365 / M1 0.2332. No lactate/regulator gain.
- No formal outputs overwritten; no raw big files added.
