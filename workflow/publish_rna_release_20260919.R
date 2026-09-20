#!/usr/bin/env Rscript
# Register a validated, visually reviewed release; no desktop writes.
suppressPackageStartupMessages({library(data.table); library(digest)})
args<-commandArgs(TRUE);stopifnot(length(args)==5L)
roles<-c('expression','qsmooth','panel','assisted','delivery')
release<-data.table(Role=roles,Path=args)
delivery<-args[[5]]
v<-fread(file.path(delivery,'validation.csv'))
visual<-fread('audit/20260919_rna_repair/visual_qa.csv')
stopifnot(all(v$Pass),nrow(visual)==18L,all(visual$Pass))
sm<-fread(file.path(args[[2]],'summary.csv'),header=FALSE,col.names=c('Metric','Value'))
fwrite(sm,file.path(delivery,'normalization_summary.csv'))
file.copy('audit/20260919_rna_repair/visual_qa.csv',file.path(delivery,'visual_qa.csv'),overwrite=TRUE)
file.copy(c('docs/RNA_DATA_PROCESSING_AND_ANALYSIS_WORKFLOW.md','docs/RNA_WORKFLOW_BRIEF.md'),delivery,overwrite=TRUE)
writeLines(c('# RNA release 2026-09-19','',
 'Current source conversion corrected; 28 distinct RNA references / 1,898 unique sample IDs / 31 proteome mappings.',
 'Statistical models are unchanged at the user\'s request. Updated inputs are evaluated by the existing methods.',
 '', '## Deliverables',
 '* `core/`: 14 PDF/PNG pairs plus underlying tables; frozen 196-gene G2M panel, fixed 371-gene DDR panel, 5,224 Kla genes, 48 regulators (49 role entries).',
 '* `reference/`: 4 PDF/PNG pairs plus three-source descriptive tables; 31 rows retain explicit shared-reference mapping.',
 '* `source_changes.csv`: only BPH and TALL-104 source expression values change; all other 26 references remain numerically unchanged.',
 '* `validation.csv`, `visual_qa.csv`: numerical/scope and rendered-output checks.',
 '* `release_sha256.csv`: manifest of this release and its supporting matrices, metadata and scripts (paths relative to repository root).',
 '* `normalization_summary.csv`: separate mean and maximum A/B differences; A and B weights reported separately.',
 '', 'Historical releases are retained for audit and are not merged into this directory.',
 'Source/batch and donor/model limitations remain as documented in the Methods. Passing these checks is not statistical-method endorsement.'),
 file.path(delivery,'README.md'))
fwrite(release,'config/rna_current_release.csv')
inputs<-c(list.files('data/publication_input',full.names=TRUE,recursive=TRUE),
 'config/rnaseq_expression_extraction_contract_20260919.csv',
 'config/rna_hallmark_g2m_checkpoint_frozen_20260919.csv',
 'audit/20260916_full_31_rna_status/rna_31_group_status.csv')
inputs<-sort(inputs[file.exists(inputs)&!dir.exists(inputs)&basename(inputs)!='.DS_Store'])
hashes<-function(paths)data.table(Path=paths,SHA256=vapply(paths,function(f)digest(file=f,algo='sha256'),''))
fwrite(hashes(inputs),'config/kla31_current_inputs_sha256.csv')
files<-unlist(lapply(args,list.files,full.names=TRUE,recursive=TRUE))
files<-files[file.exists(files)&!dir.exists(files)&basename(files)!='.DS_Store'&!grepl('release_sha256.csv$',files)]
scripts<-c('build_31_expression_matrices_20260916.R','lib_kla31_expression_20260916.R',
 'finalize_31_expression_matrices_20260916.R','verify_31_expression_matrices_20260916.R',
 'qsmooth_31group_20260916.R','overlay_ddr_panel_20260916.py','explore_rna_assisted_ddr_20260917.R',
 'plot_rna_reference_31group_20260917.R','render_teacher_questions_rna_exact.R',
 'rebuild_rna_release_20260919.sh','rebuild_rna_server_20260919.sh',
 'validate_rna_release_20260919.R','publish_rna_release_20260919.R','preflight_kla31.R',
 'fingerprint_rna_sources_20260919.R','sync_rna_delivery_20260919.py',
 'tests/test_expression_conversion_20260919.R')
files<-sort(unique(c(files,file.path('workflow',scripts),'config/rna_current_release.csv',
                    'config/kla31_current_inputs_sha256.csv')))
fwrite(hashes(files),file.path(delivery,'release_sha256.csv'))
cat('RELEASE_REGISTERED',delivery,'files',length(files),'\n')
