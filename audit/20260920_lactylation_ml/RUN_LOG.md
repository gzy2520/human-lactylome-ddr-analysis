# RUN_LOG — 20260920 lactylation ML feasibility (UTC 2026-09-20)

Seed: 25. R 4.4.3 (aarch64). Packages: data.table/digest/ggplot2 (base stats for simulation).
Note: glmnet/caret/randomForest/xgboost/pROC/PRROC NOT installed locally -> also blocks formal training.

- 22:06 mkdir workflow/lactylation_ml config/lactylation_ml outputs/20260920_lactylation_ml audit/20260920_lactylation_ml
- 22:06 wrote config/lactylation_ml/task_feasibility.csv (6 tasks A1/A2/A3/B1/B2/C; B1 audit-only)
- 22:06 wrote config/lactylation_ml/feature_dictionary.csv (7 rows; lactate-axis NOT frozen)
- 22:06 wrote config/lactylation_ml/frozen_params.csv (seed 25; qsmooth descriptive-only)
- 22:08 ran 01_build_observation_contract.R -> observation_contract.csv (31 rows) + connectivity_groups.csv (initial 28-key version)
- 22:09 FIX 01b_fix_contract.R: KLA31_03/04 linkage -> same_study_different_RNA_ref; connectivity -> 19 components over shared-RefKey OR shared-PXD (sizes 1,1,2,1,1,1,2,2,1,4,5,1,1,1,2,1,1,2,1)
- 22:09 ran 02_gate_and_simulation.R -> training_gate.csv (7 FAIL) + simulation_quarantine_DO_NOT_REPORT/ (19 groups synthetic; ridge lambda=1 LOGO plumbing OK)
- 22:10 ran 03_descriptive_audit.R -> audit_detection_confounding.csv (spearman 0.189 p=0.308) + audit_by_category.csv + audit_matched_coverage.csv (21 rows/3 PXD) + audit_regulator_spread.csv (49 entries) + audit_connectivity.csv + fig_audit_detection_confounding.{png,pdf} + sourcedata
- 22:10 refreshed audit/20260920_lactylation_ml/input_hashes.csv (SHA-256 of 19 non-quarantine files)
- Preflight: workflow/preflight_kla31.R FAILS locally (missing outputs/20260919_expression_corrected/* server-only matrices). NOT bypassed; this delivery uses only local frozen inputs. Must be disclosed in commit message and REPORT §3.

Formal training: BLOCKED. No real model performance produced. Simulation numbers are synthetic and quarantined.
