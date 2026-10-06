#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7L)
sample <- args[1]; arc <- args[2]; bed <- args[3]; json <- args[4]
set.seed(123)
cfg <- read_analysis_config(json, sample)
limits <- cfg$thresholds
assays <- make_assays(arc, bed)
tab <- qc_table(sample, assays)
keep <- apply_qc(tab, limits)
if (length(keep) < 100L) stop_bad("Fewer than 100 joint QC cells: ", sample)
rna <- subset(assays$rna, cells = keep)
rna <- preprocess_qc_rna(rna, cfg$clustering)
rna <- remove_qc_clusters(rna, cfg$doublet$exclude_clusters)
cluster_pass <- ncol(rna)
# Recompute features/PCA/neighbours after removing reviewed low-quality clusters.
rna <- preprocess_qc_rna(rna, cfg$clustering)
dims <- cfg$clustering$RNA_dims
suppressPackageStartupMessages(library(DoubletFinder))
sweep <- paramSweep(rna, PCs = dims, sct = FALSE)
bcmvn <- find.pK(summarizeSweep(sweep, GT = FALSE))
pk <- suppressWarnings(as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])))
if (!is.finite(pk)) stop_bad("DoubletFinder pK sweep failed: ", sample)
n_exp <- max(1L, round(ncol(rna) * cfg$doublet$expected_rate *
  (1 - modelHomotypic(rna$pre_doublet_clusters))))
rna <- doubletFinder(rna, PCs = dims, pK = pk, nExp = n_exp, sct = FALSE)
df_col <- grep("^DF.classifications", colnames(rna[[]]), value = TRUE)
if (length(df_col) != 1L) stop_bad("Ambiguous DoubletFinder result: ", sample)
singlets <- colnames(rna)[rna[[]][[df_col]] == "Singlet"]
if (length(singlets) < 50L) stop_bad("Fewer than 50 singlets: ", sample)
rna <- subset(rna, cells = singlets)
atac <- subset(assays$atac, cells = singlets)
atac <- atac[, colnames(rna)]
rna$barcode <- colnames(rna); atac$barcode <- colnames(atac)
rna$sample_id <- sample; atac$sample_id <- sample
rna <- RenameCells(rna, add.cell.id = sample)
atac <- RenameCells(atac, add.cell.id = sample)
stopifnot(identical(colnames(rna), colnames(atac)))
saveRDS(rna, args[5]); saveRDS(atac, args[6])
write_tsv(data.frame(sample_id = sample, raw_joint = nrow(tab), qc_pass = length(keep),
  cluster_pass = cluster_pass, excluded_clusters = paste(cfg$doublet$exclude_clusters, collapse = ","),
  df_singlets = length(singlets), pK = pk, expected_rate = cfg$doublet$expected_rate,
  expected_doublets = n_exp), args[7])
