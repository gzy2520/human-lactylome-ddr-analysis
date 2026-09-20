#!/usr/bin/env Rscript
# 01b_fix_contract.R — patch linkage + connectivity after review.
suppressPackageStartupMessages({library(data.table)})
set.seed(25)
fp <- "outputs/20260920_lactylation_ml/observation_contract.csv"
obs <- fread(fp)
# Fix KLA31_03/04: same PXD + same study samples but DIFFERENT RNA refs -> same_study_different_RNA_ref
obs[GroupID %in% c("KLA31_03", "KLA31_04"), Linkage := "same_study_different_RNA_ref"]
# Rewrite connectivity as leakage-connected components over shared ReferenceKey AND shared PXD.
# Two rows are linked if they share ReferenceKey OR share PXD; components via BFS.
n <- nrow(obs)
adj <- matrix(FALSE, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) {
  if (obs$ReferenceKey[i]==obs$ReferenceKey[j] || obs$PXD[i]==obs$PXD[j]) adj[i, j] <- TRUE
}
comp <- rep(0L, n); cur <- 0L
for (i in seq_len(n)) if (comp[i]==0L) {
  cur <- cur + 1L; q <- i; comp[i] <- cur
  while (length(q)>0L) { v <- q[1L]; q <- q[-1L]
    nb <- which(adj[v, ] & comp==0L); comp[nb] <- cur; q <- c(q, nb) }
}
obs[, ConnGroup := comp]
conn <- obs[, .(n_rows=.N, members=paste(sort(GroupID), collapse=";"),
  refs=paste(sort(unique(ReferenceKey)), collapse=";"),
  pxds=paste(sort(unique(PXD)), collapse=";")), by=ConnGroup][order(ConnGroup)]
fwrite(obs, fp)
fwrite(conn, "outputs/20260920_lactylation_ml/connectivity_groups.csv")
fwrite(obs[, .(GroupID, PXD, SampleGroup, ReferenceKey, Linkage, SharedReference, ConnGroup)],
  "outputs/20260920_lactylation_ml/observation_contract_min.csv")
cat(sprintf("FIX_OK: linkage fixed; connected components=%d (rows: %s)\n",
  nrow(conn), paste(conn$n_rows, collapse=",")))
print(conn)
