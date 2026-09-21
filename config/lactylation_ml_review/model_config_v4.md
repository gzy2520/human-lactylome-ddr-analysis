# B-R v4: correctness repair and biological module correction
Frozen before v4 training. Seed 25. Same 6744 observations, labels, 19 connected groups, models and ridge grid as v3. No selection based on new performance.

- Inner selection: calculate validation LogLoss for EACH connected group; each training connected group contributes once to the candidate-lambda mean. This matches the outer equal-group estimand. Training likelihood remains row-weighted, unchanged from v3.
- Primary module: mean of source-material median log2(TPM+1) for four distinct L-lactate/pyruvate conversion genes: ENSG00000134333 (LDHA, Entrez 3939), ENSG00000111716 (LDHB), ENSG00000166796 (LDHC), ENSG00000171989 (LDHAL6B). This is a bidirectional expression module, NOT lactate production or flux.
- ENSG00000288299 is also mapped to Entrez 3939 / UniProt P00338 in the frozen mapping. Excluded from the primary module to avoid counting LDHA twice; retain both raw expression rows and compare them. Choose ENSG00000134333 using the current NCBI Gene 3939 reference mapping, not performance.
- LDHD (ENSG00000166816) is retained separately as D-lactate oxidation, not included in L module or fitted as a new predictor. This is a biological correction, not feature search.
- Two predeclared runs: weight_only uses the exact v3 mixed feature with repaired inner weighting; corrected_L uses the four-gene L-conversion module and repaired weighting. This separates numerical and membership changes.
- No additional algorithms/modules. Model names and solver remain as v3. Internal column mact_generation_lactate retained for API compatibility only; v4 displayed names must use L-lactate/pyruvate conversion.
- S1 and true S2 for primary corrected_L; S2 covers MP, ME_P, ME2P with the same tuned engine, all 19 deleted groups. Other models have main outer validation only.
- Report mean paired change AND actual improved/worse/tied group counts (tolerance 1e-12). Do not infer coefficient sign means uniform predictive benefit/harm.
- Primary LogLoss, AP/Brier/ROC and material calibration secondary. Same-data exploratory repair, not independent confirmation. Isolation tests describe tested cases, not a universal proof of zero leakage.

Sources checked 2026-09-21:
https://www.ncbi.nlm.nih.gov/gene/3939
https://www.uniprot.org/uniprotkb/Q9BYZ2/entry (LDHAL6B annotation; not equivalent to direct activity measurement here)
https://pubmed.ncbi.nlm.nih.gov/37863926/ (LDHD D-lactate oxidation)
https://pubmed.ncbi.nlm.nih.gov/30931947/ (human LDHD deficiency)
Frozen supporting mapping: config/lactate_metabolism/go_current_uniprot_mapping.csv.
