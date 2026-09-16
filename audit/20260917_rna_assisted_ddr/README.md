# RNA-assisted view of the DDR proteome (exploratory, 2026-09-17)

The proteome records only which proteins were captured, so an absent protein is ambiguous.
The 31-group qsmooth profile supplies a continuous value per gene and separates
'not expressed' from 'expressed but not captured as lactylated'.

Exploratory: the RNA reference is a material-class profile drawn from different studies than
the proteome, so groups are compared as classes, not as paired samples.
`expression_threshold_sweep.csv` and `within_group_percentile_sweep.csv` exist because the
headline share depends on where the expression cut is placed.
