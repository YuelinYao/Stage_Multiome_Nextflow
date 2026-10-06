#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 8L)
sample <- args[1]; arc <- args[2]; bed <- args[3]; json <- args[4]
set.seed(123)
cfg <- read_pipeline_config(json, sample)
arc_detected_cells <- read_arc_detected_cells(arc)
rate_source <- if (identical(cfg$doublet$expected_rate, "auto"))
  "ARC summary; Next GEM Multiome table" else "Pipeline JSON override"
if (identical(cfg$doublet$expected_rate, "auto"))
  cfg$doublet$expected_rate <- nextgem_multiome_doublet_rate(arc_detected_cells)
message(sample, ": ARC detected cells = ", arc_detected_cells,
        "; expected doublet rate = ", cfg$doublet$expected_rate, " (", rate_source, ")")
limits <- cfg$thresholds
assays <- make_assays(arc, bed)
tab <- qc_table(sample, assays)
keep <- apply_qc(tab, limits)
message(sample, ": initial QC retained ", length(keep), " / ", nrow(tab),
        " cells (removed ", nrow(tab) - length(keep), ")")
if (length(keep) < 100L) stop_bad("Fewer than 100 joint QC cells: ", sample)
rna <- subset(assays$rna, cells = keep)
# Use the saved barcode-to-cluster assignments from 05a. No preliminary
# clustering is repeated here: exclusions refer to the exact reviewed labels.
reviewed_clusters <- read_reviewed_qc_clusters(args[8], sample, tab, colnames(rna))
rna$reviewed_qc_cluster <- unname(reviewed_clusters[colnames(rna)])
rna$pre_doublet_clusters <- rna$reviewed_qc_cluster
rna <- remove_qc_clusters(rna, cfg$doublet$exclude_clusters)
cluster_pass <- ncol(rna)
message(sample, ": reviewed-cluster removal retained ", cluster_pass, " / ", length(keep),
        " cells (removed ", length(keep) - cluster_pass, ")")
# Recompute features/PCA/neighbours after removing reviewed low-quality clusters.
# Use standard RunPCA defaults (50 computed PCs) and the doublet.PC_dims range.
# These DoubletFinder PCs are independent of stage/review RNA_dims settings.
dims <- cfg$doublet$PC_dims
doublet_clustering <- cfg$doublet
doublet_clustering$RNA_dims <- dims
rna <- preprocess_qc_rna(rna, doublet_clustering, default_pca = TRUE)
message(sample, ": DoubletFinder uses PCs ", min(dims), "-", max(dims),
        " (", ncol(Embeddings(rna, "pca")),
        " PCs computed)")
suppressPackageStartupMessages(library(DoubletFinder))
sweep <- paramSweep(rna, PCs = dims, sct = FALSE)
bcmvn <- find.pK(summarizeSweep(sweep, GT = FALSE))
bcmvn <- bcmvn[is.finite(bcmvn$BCmetric), , drop = FALSE]
if (!nrow(bcmvn)) stop_bad("No finite DoubletFinder BCmetric: ", sample)
pk <- suppressWarnings(as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])))
if (length(pk) != 1L || !is.finite(pk) || pk <= 0 || pk >= 1) stop_bad("DoubletFinder pK sweep failed: ", sample)
homotypic_proportion <- modelHomotypic(rna$pre_doublet_clusters)
expected_doublets_unadjusted <- round(ncol(rna) * cfg$doublet$expected_rate)
n_exp <- max(1L, round(ncol(rna) * cfg$doublet$expected_rate *
  (1 - homotypic_proportion)))
rna <- doubletFinder(rna, PCs = dims, pK = pk, nExp = n_exp, sct = FALSE)
df_col <- grep("^DF.classifications", colnames(rna[[]]), value = TRUE)
if (length(df_col) != 1L) stop_bad("Ambiguous DoubletFinder result: ", sample)
classifications <- as.character(rna[[]][[df_col]])
if (anyNA(classifications) || any(!classifications %in% c("Singlet", "Doublet")))
  stop_bad("Invalid DoubletFinder classifications: ", sample)
# Plot all cells assessed by DoubletFinder before selecting singlets. UMAP
# uses the same configured PCs as the sweep and doublet classification.
rna$DoubletFinder <- factor(classifications, levels = c("Singlet", "Doublet"))
rna <- RunUMAP(rna, reduction = "pca", dims = dims, seed.use = 123, verbose = FALSE)
p <- DimPlot(rna, reduction = "umap", group.by = "DoubletFinder") +
  ggplot2::ggtitle(paste0(sample, ": DoubletFinder")) +
  ggplot2::theme_bw() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
ggplot2::ggsave(filename = file.path(dirname(args[7]), paste0(sample, "_Doublet.pdf")),
  plot = p, width = 6, height = 5, units = "in")
singlets <- colnames(rna)[classifications == "Singlet"]
message(sample, ": DoubletFinder retained ", length(singlets), " / ", cluster_pass,
        " singlets (removed ", cluster_pass - length(singlets), ")")
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
  arc_detected_cells = arc_detected_cells, rate_source = rate_source,
  qc_removed = nrow(tab) - length(keep),
  cluster_pass = cluster_pass, excluded_clusters = paste(cfg$doublet$exclude_clusters, collapse = ","),
  cluster_removed = length(keep) - cluster_pass,
  df_singlets = length(singlets), df_removed = cluster_pass - length(singlets),
  df_pcs = paste(dims, collapse = ","), pca_computed = ncol(Embeddings(rna, "pca")),
  pK = pk, expected_rate = cfg$doublet$expected_rate,
  homotypic_proportion = homotypic_proportion,
  expected_doublets_unadjusted = expected_doublets_unadjusted,
  expected_doublets = n_exp), args[7])
