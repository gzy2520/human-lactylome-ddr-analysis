# Integrated human lysine lactylation proteomics and DNA-damage response

This repository contains the publication workflow for a cross-study human
lysine lactylation (Kla) proteomics analysis. The final comparison uses 33
sample groups with whole-proteome references and a fixed set of 507
Kla-intersecting DDR proteins.

## Reproduce

```bash
python3 -m pip install -r requirements.txt
Rscript workflow/install_r_dependencies.R
Rscript workflow/run_pipeline.R all
Rscript workflow/record_environment.R .
Rscript workflow/build_manifest.R .
```

Raw and processed PXD files are kept under `data/` locally and are not tracked
by Git. Analysis code uses isoform-stripped UniProt `BaseAccession`; gene
symbols are display and audit annotations only. All stochastic analyses use
seed 25. Exact R versions used for each run are written to
`results/provenance/`; `DESCRIPTION` declares the required packages.

See `PROJECT_INDEX.md` for the directory map and
`NEW_CHAT_PROJECT_PROMPT.md` for the full analysis contract.
