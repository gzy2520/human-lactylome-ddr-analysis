#!/usr/bin/env Rscript
p <- 'audit/20260915_replicate_resolution/resolved_matrices'
inv <- rbind(read.csv(file.path(p,'resolved_reference_inventory.csv')),read.csv(file.path(p,'native_reference_inventory.csv')))
samples <- rbind(read.csv(file.path(p,'selected_samples.csv')),read.csv(file.path(p,'native_selected_samples.csv'))[,c('ReferenceKey','GSE','GSM','Condition')])
gate <- read.csv('config/rnaseq_teacher_review_gate_20260915.csv')
gate$MatrixSampleN <- NA_integer_; gate$NumericQCPass <- FALSE
gate$StableID <- ''; gate$ConditionCaveat <- ''; gate$AnalysisReady <- FALSE
for(i in seq_len(nrow(inv))) {
 s <- inv[i,]; idx <- which(gate$GroupID %in% strsplit(s$GroupIDs,';',fixed=TRUE)[[1]])
 gate$ProposedReplacementSource[idx] <- s$GSE
 gate$ProposedReplacementSamples[idx] <- paste(samples$GSM[samples$ReferenceKey==s$ReferenceKey],collapse=';')
 gate$ProposedReplacementStatus[idx] <- 'stable-ID count matrix validated; biological eligibility review remains'
 gate$MatrixSampleN[idx] <- s$MatrixN; gate$NumericQCPass[idx]<-TRUE
 gate$StableID[idx] <- s$StableID; gate$ConditionCaveat[idx] <- s$Condition
 gate$ReviewAction[idx] <- 'use new candidate for source-local QC; retain teacher color and condition caveat'
}
stopifnot(nrow(gate)==31,!anyDuplicated(gate$GroupID),sum(gate$NumericQCPass)==13)
write.csv(inv,file.path(p,'combined_reference_inventory.csv'),row.names=FALSE)
write.csv(gate,'config/rnaseq_replicate_resolution_gate_20260915.csv',row.names=FALSE)
