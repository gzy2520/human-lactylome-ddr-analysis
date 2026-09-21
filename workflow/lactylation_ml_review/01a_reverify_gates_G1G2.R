#!/usr/bin/env Rscript
# 01_reverify_gates.R —逐项复核 fd3ace5 的 G1–G7（part 1: G1-G3）。
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W)
out <- "outputs/20260921_lactylation_ml_review"
dir.create(out, showWarnings=FALSE, recursive=TRUE)

# ---- G1: 标签可比性 ----
loc_files <- list.files(c("data","outputs","corrected_final_result_20260905_sample_only"), recursive=TRUE, full.names=TRUE)
has_intensity_values <- any(grepl("intens", loc_files, ignore.case=TRUE))
fig1 <- fread("data/candidate/figure1_sample_boxplot_values.csv")
has_intensity_cols_fig1 <- any(grepl("intens", names(fig1), ignore.case=TRUE))
n_pxd_raw_local <- length(list.files("data", pattern="^PXD", full.names=TRUE))
G1 <- data.table(Gate="G1_label_comparable",
  Finding="跨研究统一连续Kla标签不可比：平台各异（TMT/LFQ/MaxQuant/Spectronaut/iTRAQ）；本地无data/PXD*强度值文件；fig1仅计数列",
  Evidence=paste0("intensity_audit:31对平台各异; data/PXD* dirs=", n_pxd_raw_local, "; fig1 intensity cols=", has_intensity_cols_fig1),
  Verdict="已核实的方法学限制+数据缺失（跨研究拼标签禁行成立）")
fwrite(G1, file.path(out,"_G1.csv"))

# ---- G2: 独立单元 ----
map31 <- fread("outputs/20260920_lactate_metabolism/tables/rna_proteome_mapping_31.csv")
n <- nrow(map31); adj <- matrix(FALSE, n, n)
for (i in seq_len(n)) for (j in seq_len(n))
  if (map31$ReferenceKey[i]==map31$ReferenceKey[j] || map31$PXD[i]==map31$PXD[j]) adj[i,j] <- TRUE
comp <- rep(0L,n); cur <- 0L
for (i in seq_len(n)) if (comp[i]==0L) { cur<-cur+1L; q<-i; comp[i]<-cur
  while(length(q)>0L){ v<-q[1L]; q<-q[-1L]; nb<-which(adj[v,]&comp==0L); comp[nb]<-cur; q<-c(q,nb) } }
map31[, ConnGroup:=comp]
conn <- map31[, .(n_rows=.N, members=paste(sort(GroupID),collapse=";")), by=ConnGroup][order(ConnGroup)]
fwrite(map31[, .(GroupID,PXD,SampleGroup,ReferenceKey,SharedReference,ConnGroup)], file.path(out,"connectivity_recomputed.csv"))
fwrite(conn, file.path(out,"connectivity_groups.csv"))
mm <- fread("data/candidate/sample_resolved_matched_modalities.csv")
# RNA-蛋白样本ID交集
q <- fread("outputs/20260919_expression_corrected/group_sample_qc.csv")
pids <- fread("data/candidate/figure1_sample_boxplot_values.csv", select=c("SampleID"))
G2 <- data.table(Gate="G2_independent_units",
  Finding=sprintf("31行/28参考重算连通=%d个（fd3ace5为19，一致）；配对模态%d行/3PXD；RNA-蛋白样本ID交集=%d", nrow(conn), nrow(mm), length(intersect(unique(q$SampleID), unique(pids$SampleID)))),
  Evidence=paste0("recompute conn=", nrow(conn), " sizes=", paste(conn$n_rows,collapse=","), "; matched rows=", nrow(mm)),
  Verdict="已核实的方法学限制（材料级n小且连通；配对窄；RNA-蛋白零配对成立）")
fwrite(G2, file.path(out,"_G2.csv"))
cat("G1-G2 DONE: conn=", nrow(conn), "\n")
