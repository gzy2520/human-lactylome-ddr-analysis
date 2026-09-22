# B-R v5 execution log (2026-09-22)

- Started from clean a08f115 on research/kla31-rnaseq-20260912. Existing RNA preflight passed (163 frozen / 17 current inputs / 174 release files); no RNA inputs or released results changed.
- Inspected v4 group errors; specified nine candidate frequency configurations and two adaptive models in v5_contract.md before evaluating their outer performance. No expansion after viewing results. Hashes were recorded after execution; they are integrity records, not externally timestamped preregistration.
- R engine tests passed before training. run.R completed all 19 outer folds; evaluation preserved v4 fixed comparisons. Original MP reproduced to <1e-12. Group-weighted calibrated probabilities improved mean LogLoss slightly but worsened material calibration MAE.
- Added post-result robustness check without model changes: sensitivity.R deleted each of 19 groups and repeated full nested training for the remaining 18 outer folds (342 fits). Kept contemporaneously recomputed MP as the S2 comparator; did not claim this reproduces S2 MP_cal.
- Post-training tests.R and verify_results.R: 27 checks PASS. Independent base-R group-loss reconstruction matches results. 103 historical v4 manifest files remain unchanged.
- evaluate.R rendered two PDF/PNG pairs; inspected both PDFs via pdftoppm. No clipping or overlapping labels; calibration plots use identical [0,1] axes. Paired-loss plots have different y scales, visible on axes, to show each model's fold differences.
- Fixed a syntax-only extra parenthesis in report.R; no training or metric changes. Generated Chinese report directly from result tables.
- Reproduction: use an isolated copy, remove only that copy's v5 outputs, then bash workflow/lactylation_ml_optimize/run_all.sh. Hash packaging and PDF visual QA are separate steps.
- Input and deliverable SHA-256 manifests accompany this log. Server delivery, when completed, is documented in server_verification.log; manifest excludes itself and that verification log to avoid self-reference.
