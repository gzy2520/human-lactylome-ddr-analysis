#!/usr/bin/env Rscript
# Descriptive cross-cohort bridge. Stable IDs are join keys; symbols are labels only.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
})
set.seed(25)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: analyze.R OUTPUT_DIR RAW_REPOSITORY_ROOT")
out <- args[[1L]]
raw_root <- args[[2L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop("Refusing nonempty output directory")
}
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
put <- function(x, name) fwrite(x, file.path(out, "tables", paste0(name, ".csv")))

input <- c(
  script = "workflow/rna_protein_bridge/analyze.R",
  rna = "outputs/20260920_lactate_metabolism/tables/focus_expression_differences.csv",
  mapping = "config/lactate_metabolism/curated_mechanism_evidence.csv",
  group = "data/publication_input/group_summary_31.csv",
  kla = "data/publication_input/kla_protein_membership_31.csv",
  protein = file.path(raw_root, "data/PXD066054/search_results/extracted/DA/Protein_Quant.tsv"),
  site = file.path(raw_root, "data/PXD066054/search_results/extracted/PLa/XB08700B1DPLa_-PTMSiteReport.tsv")
)
if (!all(file.exists(input))) stop("Missing input: ", paste(input[!file.exists(input)], collapse = ", "))
manifest <- data.table(Role = names(input), Path = normalizePath(input),
                       SHA256 = vapply(input, digest::digest, character(1),
                                       algo = "sha256", file = TRUE, USE.NAMES = FALSE))
put(manifest, "input_manifest")

# The 21-gene focus panel was fixed in the first-stage lactate analysis.
focus <- c("LDHA", "LDHB", "SLC16A1", "SLC16A3", "MPC1", "MPC2", "PDHA1", "PDK1",
           "PC", "PCK1", "PCK2", "FBP1", "G6PC1", "ACSS2", "SUCLG1", "SUCLG2",
           "AARS1", "AARS2", "GLO1", "HAGH", "LDHD")
mapping <- unique(fread(input[["mapping"]])[Symbol %chin% focus & Measured == TRUE,
                                           .(Ensembl, Symbol, UniProt)])
stopifnot(nrow(mapping) == length(focus), !anyDuplicated(mapping$Ensembl),
          !anyDuplicated(mapping$UniProt), setequal(mapping$Symbol, focus))
mapping[, Order := match(Symbol, focus)]
rna <- fread(input[["rna"]])[Ensembl %chin% mapping$Ensembl,
                              .(Ensembl, Comparison, RNAReference = ReferenceA,
                                RNATumor = TumorA, RNADelta = DeltaA)]
stopifnot(nrow(rna) == 42L, !anyDuplicated(rna, by = c("Ensembl", "Comparison")))

groups <- fread(input[["group"]])
kla <- fread(input[["kla"]])
pair <- data.table(
  Comparison = c("HCC vs adjacent liver", "Prostate cancer vs BPH"),
  PXD = c("PXD075377", "PXD066054"),
  ReferenceGroup = c("adjacent liver", "BPH"),
  TumorGroup = c("HCC", "prostate cancer")
)
# Verify the exact sample-group strings against the frozen group summary.
for (i in seq_len(nrow(pair))) {
  actual <- groups[PXD == pair$PXD[i], SampleGroup]
  if (!all(c(pair$ReferenceGroup[i], pair$TumorGroup[i]) %in% actual)) {
    stop("Group key mismatch for ", pair$Comparison[i], ": ", paste(actual, collapse = ", "))
  }
}
bridge <- merge(mapping, rna, by = "Ensembl", allow.cartesian = FALSE)
for (i in seq_len(nrow(pair))) {
  idx <- bridge$Comparison == pair$Comparison[i]
  ref <- unique(kla[PXD == pair$PXD[i] & SampleGroup == pair$ReferenceGroup[i], BaseAccession])
  tum <- unique(kla[PXD == pair$PXD[i] & SampleGroup == pair$TumorGroup[i], BaseAccession])
  bridge[idx, `:=`(KlaReferenceDetected = UniProt %chin% ref,
                   KlaTumorDetected = UniProt %chin% tum)]
}
stopifnot(!anyNA(bridge$KlaReferenceDetected), !anyNA(bridge$KlaTumorDetected))

# PXD066054 ordinary proteome: accept only an unambiguous, exact UniProt protein group.
# Multimapped PG.ProteinGroups (semicolon-separated) are reported as missing, not assigned.
protein <- fread(input[["protein"]])
stopifnot(!anyDuplicated(protein$PG.ProteinGroups))
qcols <- grep("\\.PG\\.Quantity$", names(protein), value = TRUE)
ref_cols <- qcols[grepl("_NAT[1-5]_", qcols)]
tum_cols <- qcols[grepl("_PCa[1-5]_", qcols)]
stopifnot(length(qcols) == 10L, length(ref_cols) == 5L, length(tum_cols) == 5L)
pq <- merge(mapping, protein[, c("PG.ProteinGroups", qcols), with = FALSE],
            by.x = "UniProt", by.y = "PG.ProteinGroups", all.x = TRUE, sort = FALSE)
stopifnot(nrow(pq) == length(focus))
pq[, `:=`(ProteinExactGroup = rowSums(!is.na(as.matrix(.SD))) > 0L), .SDcols = qcols]
pq[, `:=`(ProteinBPHN = rowSums(is.finite(as.matrix(.SD)) & as.matrix(.SD) > 0)), .SDcols = ref_cols]
pq[, `:=`(ProteinPCaN = rowSums(is.finite(as.matrix(.SD)) & as.matrix(.SD) > 0)), .SDcols = tum_cols]
pq[, ProteinBPHMedian := apply(.SD, 1, function(v) if (sum(is.finite(v) & v > 0) >= 3L) median(v[is.finite(v) & v > 0]) else NA_real_), .SDcols = ref_cols]
pq[, ProteinPCaMedian := apply(.SD, 1, function(v) if (sum(is.finite(v) & v > 0) >= 3L) median(v[is.finite(v) & v > 0]) else NA_real_), .SDcols = tum_cols]
pq[, ProteinLog2Ratio := log2(ProteinPCaMedian / ProteinBPHMedian)]
protein_summary <- pq[, .(Ensembl, UniProt, ProteinExactGroup, ProteinBPHN, ProteinPCaN,
                          ProteinBPHMedian, ProteinPCaMedian, ProteinLog2Ratio)]
bridge <- merge(bridge, protein_summary, by = c("Ensembl", "UniProt"), all.x = TRUE, sort = FALSE)
bridge[Comparison == "HCC vs adjacent liver", (names(protein_summary)[-(1:2)]) := NA]
setorder(bridge, Comparison, Order)
put(bridge[, !"Order"], "fixed21_cross_layer")
put(bridge[Comparison == "Prostate cancer vs BPH", !"Order"], "prostate_cross_layer")

# Protein set overlap checks show that shared detection is a property of the source,
# while quantities and site signals may still vary.
pair_summary <- rbindlist(lapply(seq_len(nrow(pair)), function(i) {
  p <- pair[i]
  ref <- groups[PXD == p$PXD & SampleGroup == p$ReferenceGroup]
  tum <- groups[PXD == p$PXD & SampleGroup == p$TumorGroup]
  a <- unique(kla[PXD == p$PXD & SampleGroup == p$ReferenceGroup, BaseAccession])
  b <- unique(kla[PXD == p$PXD & SampleGroup == p$TumorGroup, BaseAccession])
  stopifnot(length(a) == ref$KlaProteinCount, length(b) == tum$KlaProteinCount)
  data.table(Comparison = p$Comparison, ReferenceGroup = p$ReferenceGroup, TumorGroup = p$TumorGroup,
             ReferenceKlaProteinCount = length(a), TumorKlaProteinCount = length(b),
             SharedKlaProteinCount = length(intersect(a, b)),
             KlaProteinJaccard = length(intersect(a, b)) / length(union(a, b)),
             ReferenceKlaDDRPercent = 100 * ref$KlaDdrFraction,
             TumorKlaDDRPercent = 100 * tum$KlaDdrFraction,
             ReferenceOrdinaryProteinCount = ref$ReferenceProteinCount,
             TumorOrdinaryProteinCount = tum$ReferenceProteinCount,
             FocusKlaChanged = bridge[Comparison == p$Comparison,
                                      sum(KlaReferenceDetected != KlaTumorDetected)])
}))
put(pair_summary, "pair_detection_context")

# Illustrative LDHA/LDHB site quantities. Site signals are not occupancy because
# the modification and ordinary-proteome arms are separately measured.
site <- fread(input[["site"]], select = c("R.Condition", "PTM.ProteinId", "PTM.ModificationTitle",
                                          "PTM.SiteAA", "PTM.SiteProbability", "PTM.SiteLocation",
                                          "PTM.Quantity"))
# Reproduce both frozen PXD066054 protein-union sets directly from the source.
raw_ids <- function(x) {
  ids <- unlist(strsplit(x[!is.na(x)], "[;,]"), use.names = FALSE)
  ids <- trimws(ids)
  ids <- sub("^.*\\|([^|]+)\\|.*$", "\\1", ids)
  ids <- sub("-[0-9]+$", "", ids)
  sort(unique(ids[grepl("^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]([A-Z0-9]{3}[0-9])?)$", ids)]))
}
valid_site <- site[PTM.ModificationTitle == "L-Lac(K)" & PTM.SiteAA == "K" &
                     PTM.SiteProbability > 0]
raw_bph <- raw_ids(valid_site[grepl("^NAT", R.Condition), PTM.ProteinId])
raw_pca <- raw_ids(valid_site[grepl("^PCa", R.Condition), PTM.ProteinId])
frozen_bph <- sort(unique(kla[PXD == "PXD066054" & SampleGroup == "BPH", BaseAccession]))
frozen_pca <- sort(unique(kla[PXD == "PXD066054" & SampleGroup == "prostate cancer", BaseAccession]))
stopifnot(identical(raw_bph, frozen_bph), identical(raw_pca, frozen_pca))
put(data.table(Group = c("BPH", "Prostate cancer"), RawProteins = c(length(raw_bph), length(raw_pca)),
               FrozenProteins = c(length(frozen_bph), length(frozen_pca)), ExactSetMatch = TRUE),
    "raw_kla_set_reproduction")
ldh_ids <- mapping[Symbol %chin% c("LDHA", "LDHB"), UniProt]
site <- site[PTM.ProteinId %chin% ldh_ids & PTM.ModificationTitle == "L-Lac(K)" &
             PTM.SiteAA == "K" & PTM.SiteProbability >= 0.75 & PTM.Quantity > 0]
site[, Group := fifelse(grepl("^NAT[1-5]$", R.Condition), "BPH",
                        fifelse(grepl("^PCa[1-5]$", R.Condition), "PCa", NA_character_))]
stopifnot(!anyNA(site$Group))
sample_sites <- site[, .(SiteQuantity = median(PTM.Quantity)),
                     by = .(PTM.ProteinId, PTM.SiteLocation, Group, R.Condition)]
site_summary <- sample_sites[, .(N = .N, MedianSiteQuantity = median(SiteQuantity)),
                             by = .(PTM.ProteinId, PTM.SiteLocation, Group)]
site_summary <- dcast(site_summary, PTM.ProteinId + PTM.SiteLocation ~ Group,
                      value.var = c("N", "MedianSiteQuantity"))
site_summary[, SiteLog2Ratio := log2(MedianSiteQuantity_PCa / MedianSiteQuantity_BPH)]
site_summary <- merge(site_summary, mapping[, .(UniProt, Ensembl, Symbol)],
                      by.x = "PTM.ProteinId", by.y = "UniProt")
setorder(site_summary, Symbol, PTM.SiteLocation)
put(site_summary, "ldh_site_quantities_exploratory")
sp <- copy(site_summary[N_BPH >= 3L & N_PCa >= 3L])
sp[, SiteLabel := paste0(Symbol, " K", PTM.SiteLocation)]
sp[, SiteLabel := factor(SiteLabel, levels = rev(SiteLabel))]
ps <- ggplot(sp, aes(SiteLog2Ratio, SiteLabel)) +
  geom_vline(xintercept = 0, colour = "grey75") +
  geom_segment(aes(x = 0, xend = SiteLog2Ratio, yend = SiteLabel,
                   colour = SiteLog2Ratio > 0), linewidth = 0.8) +
  geom_point(aes(colour = SiteLog2Ratio > 0), size = 2.5) +
  scale_colour_manual(values = c(`FALSE` = "#2166AC", `TRUE` = "#B2182B"), guide = "none") +
  labs(title = "LDH lactylated-site quantity is site dependent",
       subtitle = "PXD066054; median exported modified-site quantities, 3–5 detected samples per group",
       x = "Log2 ratio of site-quantity medians (prostate cancer / BPH)", y = NULL,
       caption = "Exploratory modified-site signal only; separate ordinary-proteome data do not establish site occupancy.") +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "figures", "LDH_site_quantity.pdf"), ps, width = 8, height = 5.5,
       device = cairo_pdf)
ggsave(file.path(out, "figures", "LDH_site_quantity.png"), ps, width = 8, height = 5.5,
       device = agg_png, dpi = 180)

plot_data <- rbindlist(list(
  bridge[Comparison == "HCC vs adjacent liver",
         .(Symbol, Order, Layer = "Liver RNA", Effect = RNADelta)],
  bridge[Comparison == "Prostate cancer vs BPH",
         .(Symbol, Order, Layer = "Prostate RNA", Effect = RNADelta)],
  bridge[Comparison == "Prostate cancer vs BPH",
         .(Symbol, Order, Layer = "Prostate protein", Effect = ProteinLog2Ratio)]
))
plot_data[, `:=`(Symbol = factor(Symbol, levels = rev(focus)),
                 Layer = factor(Layer, levels = c("Liver RNA", "Prostate RNA", "Prostate protein")))]
p <- ggplot(plot_data[!is.na(Effect)], aes(x = Effect, y = Symbol)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = Effect, yend = Symbol, colour = Effect > 0), linewidth = 0.75) +
  geom_point(aes(colour = Effect > 0), size = 2) +
  facet_grid(. ~ Layer) +
  scale_colour_manual(values = c(`FALSE` = "#2166AC", `TRUE` = "#B2182B"), guide = "none") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2), limits = c(-5, 3)) +
  labs(title = "Lactate-route transcripts and ordinary proteins",
       subtitle = "Fixed 21-gene panel | external RNA cohorts; prostate proteins: PXD066054, 5 BPH + 5 cancer samples",
       x = "Tumour minus reference (see metric definitions below)", y = NULL,
       caption = "RNA: difference on log2(TPM + 0.5) scale. Protein: log2 ratio of group medians of exported quantities.\nBlank protein rows lack an unambiguous exact UniProt group. No matched RNA, Kla occupancy or flux is inferred.") +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 8), panel.spacing.x = unit(0.7, "lines"))
ggsave(file.path(out, "figures", "RNA_protein_bridge.pdf"), p, width = 11, height = 8.5,
       device = cairo_pdf)
ggsave(file.path(out, "figures", "RNA_protein_bridge.png"), p, width = 11, height = 8.5,
       device = agg_png, dpi = 180)

stopifnot(nrow(bridge) == 42L, nrow(pair_summary) == 2L,
          bridge[Comparison == "Prostate cancer vs BPH", sum(ProteinExactGroup)] == 18L,
          pair_summary[Comparison == "Prostate cancer vs BPH", SharedKlaProteinCount] == 2116L,
          all(pair_summary$FocusKlaChanged == 0L))
release_files <- list.files(out, recursive = TRUE, full.names = TRUE)
release_files <- release_files[file.info(release_files)$isdir == FALSE]
release <- data.table(File = substring(release_files, nchar(out) + 2L),
                      SHA256 = vapply(release_files, digest::digest, character(1),
                                      algo = "sha256", file = TRUE, USE.NAMES = FALSE))
fwrite(release, file.path(out, "release_sha256.csv"))
cat("PASS: 42 gene-pair records, 18 exact prostate protein groups, raw set overlap and figure.\n")
