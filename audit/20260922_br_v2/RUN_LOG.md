# RUN_LOG — 20260922 B-R v2
Worktree: /Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912
Base commit: d11220e. Branch: research/kla31-rnaseq-20260912. Seed 25. R + data.table/ggplot2 (no new deps).
Nature: exploratory model improvement on the same dataset as d11220e; CV results inspected during development.

Steps (all run with Rscript --vanilla):
1. 10_metrics_tests.R — AP threshold-aggregation tests ALL PASS (constant=tie-handled; tie-order invariance; perfect ranking; single-class NA+note). PRROC not installed → independent-impl check replaced by exact hand example + reference numbers from requester (not hardcoded).
2. 06_build_BR_v2_dataset.R — contract 6744/2708 OK; 0 target genes missing; v2 set = 6744 rows.
3. 07_build_features_v2.R — mat_med 396 genes x 28 refs; NA target rows 0; 6 frozen modules coverage 1.00.
4. 08_train_BR_v2.R — 19 outer folds, all 4 glm models converged, no separation, no silent fallback; MP coverage saved.
   Equal-weight: M0 AP .3216/LL .6826; MP .6583/.5757; ME .4313/.6843; ME2 .4091/.6949; ME2P .6134/.6026; XEN .3787/.7147.
   z_mod coefficient negative in 19/19 folds (sd .0585) — module adds no predictive gain.
5. 09_sensitivity_v2.R — S1 (eval summary sensitivity): MP eqw_AP .630–.690. S2 (train-impact, re-run pipeline minus one conn, compared on remaining 18): MP .678–.738, ME2P .691–.746. Material-level |calib err|: M0 .0004, MP .0137, ME2P .0029.
6. 11_figures_v2.R — fig_v2_perfold_AP.{png,pdf}, fig_v2_material_calibration.{png,pdf}; render check OK.
7. Hashes computed over new files (hashes.csv); commit to research branch.

Key correction: d11220e M0 pooled "PR-AUC" 0.3108 was tie-split artifact → corrected AP 0.3216 (pooled). Per-fold M0 ROC = 0.500; pooled 0.3205 ROC is cross-fold baseline artifact (documented).
Frozen before looking at new model results: model_config_v2.md + modules_frozen_v2.csv (module members, directions unsigned, comparison list).
