#!/usr/bin/env Rscript
# Create a dated, non-analysis-ready RNA source-routing candidate for the 31 Kla groups.
# No expression values are transformed, merged, or statistically analysed here.
set.seed(25)

candidate_path <- 'outputs/20260914_reference_material_audit/rnaseq_reference_candidate_31.csv'
route_path <- 'outputs/20260914_reference_material_audit/rnaseq_source_routing_candidate_31.csv'

catalog <- read.csv(candidate_path, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(
  identical(sort(catalog$GroupID), sprintf('KLA31_%02d', 1:31)),
  !anyDuplicated(catalog$GroupID),
  all(catalog$AnalysisReady %in% FALSE)
)

route <- data.frame(
  GroupID=sprintf('KLA31_%02d', 1:31),
  RouteTier=c(
    'strict_raw_reprocess', 'source_release_quantification', 'source_native_gene_matrix', 'source_native_gene_matrix',
    'source_release_quantification', 'uniform_recount3_gene_counts', 'uniform_recount3_gene_counts', 'source_native_gene_matrix',
    rep('source_release_quantification', 4), rep('source_release_quantification', 4), 'source_native_gene_matrix',
    rep('source_release_quantification', 4), 'strict_raw_reprocess', 'source_native_gene_matrix', 'source_release_quantification',
    'source_native_gene_matrix', 'source_native_gene_matrix', 'source_native_gene_matrix', 'source_native_gene_matrix', 'uniform_recount3_gene_counts',
    'uniform_recount3_gene_counts', 'source_release_quantification'
  ),
  RouteSource=c(
    'GSE236746 selected ENA FASTQ', 'GTEx v8 gene reads/TPM', 'GSE181540 mRNA Expression Profiling.xlsx', 'GSE181540 mRNA Expression Profiling.xlsx',
    'GTEx v8 gene reads/TPM', 'recount3 SRP148556 G026', 'recount3 SRP014899 G026', 'GSE132714 per-sample RPKM archive',
    'GDC TCGA-LIHC STAR Counts', 'GDC TCGA-ESCA STAR Counts', 'GDC TCGA-PRAD STAR Counts', 'GDC TCGA-LIHC STAR Counts',
    'DepMap Public 24Q4 protein-coding expression', 'DepMap Public 24Q4 protein-coding expression', 'DepMap Public 24Q4 protein-coding expression', 'DepMap Public 24Q4 protein-coding expression',
    'GSE163787 4cell FPKM matrix', 'DepMap Public 24Q4 protein-coding expression', 'DepMap Public 24Q4 protein-coding expression', 'DepMap Public 24Q4 protein-coding expression', 'DepMap Public 24Q4 protein-coding expression',
    'GSE235595 selected ENA FASTQ', 'GSE266884 MES28 MCT1 core table', 'DepMap Public 24Q4 protein-coding expression',
    'GSE203529 WT untreated RSEM gene results', 'GSE275256 STAR readcount matrix', 'GSE269418 reported raw-count matrix', 'GSE269418 reported raw-count matrix',
    'recount3 SRP117021 G026', 'recount3 SRP161553 G026', 'ENCODE ENCSR000COZ RSEM quantifications'
  ),
  InputUnit=c(
    'FASTQ to one declared gene-count/TPM pipeline', 'gene reads and TPM as released', 'reported Ensembl count columns and FPKM', 'reported Ensembl count columns and FPKM',
    'gene reads and TPM as released', 'gene counts', 'gene counts', 'per-sample RPKM',
    'STAR gene-count fields', 'STAR gene-count fields', 'STAR gene-count fields', 'STAR gene-count fields',
    rep('TPMLogp1 as distributed', 4), 'FPKM', rep('TPMLogp1 as distributed', 4),
    'FASTQ to one declared gene-count/TPM pipeline', 'reported read_count columns', 'TPMLogp1 as distributed',
    'RSEM expected_count, TPM and FPKM', 'STAR gene readcount', 'reported raw-count matrix', 'reported raw-count matrix',
    'gene counts', 'gene counts', 'RSEM gene quantification'
  ),
  StableAnalysisID=c(
    'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID',
    rep('Ensembl gene ID', 4), rep('Entrez ID parsed from DepMap header', 4), 'Ensembl gene ID', rep('Entrez ID parsed from DepMap header', 4),
    'Ensembl gene ID', 'Entrez gene ID', 'Entrez ID parsed from DepMap header', 'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID',
    'Ensembl gene ID', 'Ensembl gene ID', 'Ensembl gene ID'
  ),
  WhyThisRoute=c(
    'The only released matrix is transcript-level ENST FPKM; selected raw reads are needed for gene-level values.',
    'Fixed standard release for normal lung.',
    'The paper documents paired hypertrophic-scar and adjacent full-thickness normal-skin sampling.',
    'The paper documents paired hypertrophic-scar and adjacent full-thickness normal-skin sampling.',
    'Fixed standard release for hippocampus.',
    'All 21 selected SRR accessions are present in the same recount3 gene-count matrix.',
    'All five selected sperm mRNA SRRs are present in the same recount3 gene-count matrix.',
    'The source archive contains all 18 selected BPH samples as versioned-Ensembl RPKM files.',
    'Adjacent liver must remain TCGA Solid Tissue Normal.',
    'ESCC must remain squamous-filtered TCGA-ESCA Primary Tumor.',
    'Primary prostate tumour remains TCGA-PRAD.',
    'HCC remains TCGA-LIHC Primary Tumor.',
    'Fixed shared cell-line reference release.', 'Fixed shared cell-line reference release.', 'Fixed shared cell-line reference release.', 'Fixed shared cell-line reference release.',
    'Exact TALL-104 matrix exists, but has one bulk library only.',
    'Fixed shared cell-line reference release.', 'Fixed shared cell-line reference release.', 'Fixed shared cell-line reference release.', 'Fixed shared cell-line reference release.',
    'No source matrix contains the selected DMSO PC-3M controls.',
    'Selected MES28 control columns are supplied in the study table.',
    'Fixed shared cell-line reference release.',
    'Three WT untreated HEK293T RSEM gene-result files are supplied in the GEO archive.',
    'Selected vehicle columns are supplied in the STAR count matrix.',
    'Selected normoxia solvent-control columns are supplied in the matrix.',
    'The same verified HK-2 control matrix is shared with KLA31_27.',
    'All 12 selected untreated MCF10A SRRs are present in the same recount3 gene-count matrix.',
    'All nine selected NSC SRRs are present in the same recount3 gene-count matrix; biological-model mismatch remains.',
    'Released replicate quantifications are already downloaded from ENCODE.'
  ),
  CrossSourceMerge='No; retain source units and use only a predeclared cross-source summary layer.',
  AnalysisReady=FALSE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  identical(route$GroupID, catalog$GroupID),
  !anyDuplicated(route$GroupID),
  all(route$StableAnalysisID != ''),
  all(route$AnalysisReady %in% FALSE)
)
write.csv(route, route_path, row.names=FALSE, na='')
message('Wrote ', route_path, ' with ', nrow(route), ' groups; all remain AnalysisReady=FALSE.')
