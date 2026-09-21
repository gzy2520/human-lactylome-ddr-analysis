# B-R v4 execution record

- Started from 3f6f811 with a clean worktree. Historical v2/v3 files retained.
- Reproduced 9/19 (not 15/19) worse v3 folds and unequal inner row/group risk.
- Independently read 28 source matrices: all 6744 target values and old module means match. Alternate LDHA ID is absent in every source; no double contribution occurred. First strict coverage check stopped on this missing ID; corrected coverage now explicitly records NA rather than calling it measured.
- Frozen four-gene L-conversion module, separate D-oxidation, old-feature weight-only control before model runs.
- Actual inner connected-group loss means independently verified against selected-lambda records.
- Completed both main 19-fold runs and primary S2 for MP/ME_P/ME2P: all 19 deletions, each with 18 outer folds and nested tuning.
- Inherited acceptance first stopped on a renamed test helper; corrected to shared metrics_row_v3 and reran. Final 17 inherited + 7 additional tests pass.
- Four PDF/PNG pairs rendered and visually reviewed. Increased plot margins/shortened axis text after visual review; numerical outputs unchanged.
- Report and worse/better/tied counts generated from actual per-group predictions. Historical input/output hashes checked unchanged; current RNA preflight passed.
- No change to frozen RNA/Kla pipeline, labels, task observations, existing significance models, or old results.
