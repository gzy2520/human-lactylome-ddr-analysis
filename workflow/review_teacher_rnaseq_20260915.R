#!/usr/bin/env Rscript
# Preserve teacher decisions separately from our scientific/replicate review.
set.seed(25)
out <- 'audit/20260915_teacher_review'
dir.create(file.path(out,'evidence'),recursive=TRUE,showWarnings=FALSE)
src <- read.csv('outputs/20260914_reference_material_audit/rnaseq_reference_candidate_31.csv',check.names=FALSE)
stopifnot(nrow(src)==31,!anyDuplicated(src$GroupID))
src$TeacherColor <- 'unmarked'
green <- c(1,2,5,6,7,8,22,23,25,26,29,30,31)
yellow <- c(3,4,13:21,24)
src$TeacherColor[green] <- 'green'
src$TeacherColor[yellow] <- 'yellow'
src$TeacherColor[c(27,28)] <- 'red'
src$TeacherMeaning <- c(green='passed by teacher',yellow='condition uncertain',red='incorrect per teacher',unmarked='not opened by teacher')[src$TeacherColor]
src$SelectedUnitsN <- c(3,NA,3,3,NA,21,2,18,NA,NA,NA,NA,1,1,1,1,1,1,1,1,1,3,2,1,3,3,5,5,3,9,2)
src$ObservationUnit <- 'selected GEO sample; biological independence must be confirmed'
dep <- which(src$Source=='DepMap')
src$ObservationUnit[dep] <- 'one model-level expression profile; not biological replicates'
src$ObservationUnit[c(2,5,9:12)] <- 'independent donor/case after clinical filtering and aliquot deduplication'
src$ObservationUnit[7] <- 'two donors; five runs do not increase donor n'
src$ObservationUnit[23] <- 'two MES28 shMCT1 control libraries; do not pool shKU70 controls'
src$ObservationUnit[30] <- 'nine NSC libraries; model/donor/culture hierarchy unresolved'
src$ObservationUnit[31] <- 'two released files; donor independence unverified; n cannot exceed 2'
src$ReplicateAssessment <- ifelse(is.na(src$SelectedUnitsN),'sample lock pending',ifelse(src$SelectedUnitsN<3,'below preferred n=3','nominal n>=3; verify independent units'))
src$ReviewAction <- 'retain candidate; lock sample identity/independence and complete QC'
src$ReviewAction[dep] <- 'background only; seek >=3 matched biological RNA-seq replicates; check parental versus KO/infection/co-culture'
src$ReviewAction[c(3,4)] <- 'input has no RIP antibody; exclude IP/peak tables; retain annotation/paper contradiction for review'
src$ReviewAction[c(27,28)] <- 'exclude old source from primary candidate; overnight serum starvation and mixed solvent control; seek replacement'
src$ReviewAction[c(7,17,23,31)] <- 'below preferred n=3; descriptive fallback only while seeking independent matched replicates'
src$ReviewAction[30] <- 'teacher green retained; confirm model/donor hierarchy before treating nine libraries as n=9'
src$ReviewDate <- '2026-09-15'
src$AnalysisReady <- FALSE
src$PrimaryErrorBarEligible <- FALSE
src$ProposedReplacementSource <- ''
src$ProposedReplacementSamples <- ''
src$ProposedReplacementSource[c(27,28)] <- 'GSE240748'
src$ProposedReplacementSamples[c(27,28)] <- 'GSM7708662;GSM7708663;GSM7708664'
src$ProposedReplacementStatus <- ''
src$ProposedReplacementStatus[c(27,28)] <- 'three source files inspected; gene-prefixed labels need versioned stable-ID mapping; independent culture verification pending'
write.csv(src,file.path(out,'teacher_review_31.csv'),row.names=FALSE,na='')
write.csv(src[,c('GroupID','ReferenceKey','TeacherColor','ReviewAction','PrimaryErrorBarEligible',
  'ProposedReplacementSource','ProposedReplacementSamples','ProposedReplacementStatus')],
  'config/rnaseq_teacher_review_gate_20260915.csv',row.names=FALSE)
accessions <- c('GSE240748','GSM7708662','GSM7708663','GSM7708664','GSE81194','GSM2144412','GSM2144413','GSM2144414',
 'GSE145108','GSM4494786','GSM4763919','GSE157383','GSE310291','GSM6934151','GSE65683','GSE193192')
records <- list()
for(a in accessions) {
  dest <- file.path(out,'evidence',paste0(a,'.soft.txt'))
  url <- paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',a,'&targ=self&form=text&view=full')
  code <- if(file.exists(dest) && any(startsWith(readLines(dest,warn=FALSE),'^'))) 0L else
    system2('curl',c('-fLsS','--retry','1','--max-time','20','-o',shQuote(dest),shQuote(url)))
  txt <- if(file.exists(dest))readLines(dest,warn=FALSE) else character()
  ok <- code==0 && any(startsWith(txt,paste0('^',if(startsWith(a,'GSM'))'SAMPLE' else 'SERIES',' = ',a)))
  fields <- txt[grepl('^!(Sample_(title|characteristics_ch1|treatment_protocol_ch1|growth_protocol_ch1|library_strategy|series_id|supplementary_file)|Series_(title|overall_design|sample_id|supplementary_file))',txt)]
  records[[a]] <- data.frame(Accession=a,URL=url,FetchValid=ok,Evidence=paste(fields,collapse=' || '))
}
write.csv(do.call(rbind,records),file.path(out,'new_source_evidence.csv'),row.names=FALSE)
stopifnot(sum(src$TeacherColor=='green')==13,sum(src$TeacherColor=='yellow')==12,
          sum(src$TeacherColor=='red')==2,sum(src$TeacherColor=='unmarked')==4)
cat('TEACHER_REVIEW_PASS: 31 rows; 13 green, 12 yellow, 2 red, 4 unmarked\n')
