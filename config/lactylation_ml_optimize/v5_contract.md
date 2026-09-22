# B-R v5 bounded probability optimization (2026-09-22)
Frozen before candidate performance evaluation. Seed=25. No new genes, labels, rows or algorithms beyond frequency shrinkage and two-parameter logistic probability calibration.

Task/data: same 6744 DDR-protein x material rows, 2708 positives, 19 connected groups; known-protein process detection, not biological modification truth. Primary: equal connected-group LogLoss. Secondary: AP, Brier, ROC, material calibration. All original v4 predictions retained as fixed comparison.

Motivation: v4 protein history is strongest; independent material rates vary, shared-source records receive multiple votes under row counts. Test whether shrinkage and limiting each source group's contribution improve probability estimates. No claim this resolves batch confounding.

Candidate grid (nine fixed settings):
- original: row counts, center 0.5, strength 2; unseen protein uses (training positives+1)/(training rows+2), reproducing v4 MP exactly.
- row or connected-group evidence weighting, training equal-connected-group prevalence center, strengths {0.5,2,8,32}.
- row weighting: n/k are training Ref-detected rows/positives for each protein.
- connected-group weighting: each observed protein-group contributes one vote, with mean(y) within that group as fractional positive evidence. Not a literal iid binomial posterior. s is a fixed number of prior pseudo-observations in the corresponding evidence units.
- p=(k+s*center)/(n+s). Unseen proteins use training group-mean prevalence (original as above). Clip to [1e-9,1-1e-9].

Two predeclared adaptive models:
SHRINK selects one of nine settings by inner equal-connected-group LogLoss.
SHRINK_CAL selects setting using calibrated inner validation predictions. Calibration is intercept a and slope b on logit(p), fitted with equal-connected-group row weights on cross-fitted training predictions; objective = group-mean logloss + 0.01*(b-1)^2 (intercept unpenalized), fixed penalty (not searched). No test labels in calibration.

Isolation: outer leave-connected-group-out; five inner folds. For EACH candidate/inner split, prediction on V uses only A; calibration on A uses predictions cross-fitted within A (four source-group folds). Fit final calibration using cross-fitted outer training predictions at selected setting; predict outer E with all outer train. Candidate selection only inner; all preprocessing is training-only. Store group-level losses and coefficients.

Do not use the previously viewed outer scores to revise this grid or add candidates. Include original MP and v4 MP_cal as controls. Report all candidates' adaptive results even if worse. Repeated-data exploratory development, not independent validation. Do not claim optimized to optimum, statistically significant superiority, or metabolic mechanism.

Decision: lower outer equal-group LogLoss is an observed development improvement, not a guarantee. Report paired folds and uncertainty; retain v4 as validated historical baseline. If benefit absent, stop this bounded search and report negative result.
