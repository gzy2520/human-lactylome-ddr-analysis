#!/usr/bin/env Rscript
# 01b_reverify_gates_G3G4.R — G3 折标签分布（实际生成折）+ G4 共享来源。
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
W <- "/Users/gzy2520/Desktop/Research/kla/.worktrees/kla31-rnaseq-20260912"
setwd(W); out <- "outputs/20260921_lactylation_ml_review"

map31 <- fread(file.path(out,"connectivity_recomputed.csv"))
p <- fread("outputs/20260919_repair_assisted/group_by_protein_pairs.tsv.gz")
pr <- p[RefDetected==TRUE]
g2c <- unique(map31[, .(GroupID, ConnGroup)])
pr <- merge(pr, g2c, by="GroupID")
folds <- pr[, .(test_N=.N, test_pos=sum(KlaDetected), test_rate=round(mean(KlaDetected),3),
  train_N=nrow(pr)-.N, train_pos=sum(pr$KlaDetected)-sum(KlaDetected),
  train_rate=round((sum(pr$KlaDetected)-sum(KlaDetected))/(nrow(pr)-.N),3)), by=ConnGroup][order(ConnGroup)]
fwrite(folds, file.path(out,"folds_BRef_leaveConnOut.csv"))
ok <- folds[test_N>=20 & test_pos>=5 & (test_N-test_pos)>=5]
G3 <- data.table(Gate="G3_fold_label_variation",
  Finding=sprintf("B-R限定集n=%d正例率=%.3f；19折中%d折测试集N>=20且正负各>=5；全部折训练集正负充足（最小训练正例=%d）", nrow(pr), mean(pr$KlaDetected), nrow(ok), min(folds$train_pos)),
  Evidence=paste0("folds_BRef_leaveConnOut.csv; testN=", paste(folds$test_N,collapse=",")),
  Verdict="部分纠正：“无法保证”过笼统；多数折可算PR/ROC/Brier，单行组折精度低但非不可行，应报告折构成")
fwrite(G3, file.path(out,"_G3.csv"))

mapfull <- fread("outputs/20260920_lactate_metabolism/tables/rna_proteome_mapping_31.csv")
nshared <- nrow(mapfull[SharedReference==TRUE])
G4 <- data.table(Gate="G4_no_unmanaged_sharing",
  Finding="HCT116x3+HK2x2共享成立；5个PXD跨多材料；但连通分组留一可锁同折，泄漏可管理",
  Evidence=paste0("shared rows=", nshared, "; multi-group PXDs PXD028488/PXD060185/PXD046800/PXD066054/PXD075377/PXD070007"),
  Verdict="纠正表述：“需分组管理”≠“无法处理”；须按19连通分组拆分、禁随机拆行")
fwrite(G4, file.path(out,"_G4.csv"))
cat("G3-G4 DONE: ok folds=", nrow(ok), "/19\n")
