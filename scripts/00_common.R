suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
})
stop_bad <- function(...) stop(..., call. = FALSE)
write_tsv <- function(x, file) write.table(x, file, sep = "\t", quote = FALSE,
                                           row.names = FALSE, na = "NA")
read_arc <- function(arc) {
  x <- Read10X_h5(file.path(arc, "filtered_feature_bc_matrix.h5"))
  if (!all(c("Gene Expression", "Peaks") %in% names(x)))
    stop_bad("Expected Gene Expression and Peaks in ", arc)
  x
}
annotation_mm10 <- function() {
  if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE))
    stop_bad("Install EnsDb.Mmusculus.v79 for mm10 TSS annotation")
  x <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
  GenomeInfoDb::seqlevelsStyle(x) <- "UCSC"
  GenomeInfoDb::genome(x) <- "mm10"
  x
}
make_assays <- function(arc, bed) {
  x <- read_arc(arc)
  cells <- intersect(colnames(x$`Gene Expression`), colnames(x$Peaks))
  if (!length(cells)) stop_bad("No shared RNA/ATAC barcodes: ", arc)
  rna <- CreateSeuratObject(x$`Gene Expression`[, cells, drop = FALSE],
                            min.cells = 0, min.features = 0, assay = "RNA")
  rna$percent.mt <- PercentageFeatureSet(rna, pattern = "^[Mm][Tt]-")
  rna$percent.ribo <- PercentageFeatureSet(rna, pattern = "^Rp[ls]")
  frag <- file.path(arc, "atac_fragments.tsv.gz")
  peaks <- rtracklayer::import(bed, format = "BED")
  fragments <- CreateFragmentObject(frag, cells = cells)
  mat <- FeatureMatrix(fragments = fragments, features = peaks, cells = cells,
                       keep_all_features = TRUE)
  atac <- CreateSeuratObject(CreateChromatinAssay(mat, fragments = fragments,
                              annotation = annotation_mm10(), genome = "mm10"),
                            assay = "ATAC", min.cells = 0, min.features = 0)
  atac <- NucleosomeSignal(atac, verbose = FALSE)
  atac <- TSSEnrichment(atac, verbose = FALSE)
  list(rna = rna, atac = atac)
}
qc_metric_names <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
                     "nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment")
plot_qc_violins <- function(tab, file, title) {
  if (!requireNamespace("ggpubr", quietly = TRUE))
    stop_bad("Install ggpubr to generate QC violin plots")
  ids <- unique(tab$sample_id)
  grDevices::pdf(file, width = max(10, length(ids) * 0.8), height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (metric in qc_metric_names) {
    data <- data.frame(sample_id = factor(tab$sample_id, levels = ids),
                       value = tab[[metric]])
    data <- data[is.finite(data$value), , drop = FALSE]
    if (!nrow(data)) stop_bad("No finite values for ", metric)
    p <- ggpubr::ggviolin(data, x = "sample_id", y = "value",
      fill = "sample_id", palette = "jco", trim = TRUE,
      add = "boxplot", add.params = list(width = 0.12, fill = "white")) +
      ggplot2::labs(title = paste(title, metric, sep = ": "),
                    x = "Sample", y = metric) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")
    print(p)
  }
}
qc_table <- function(sample, assays) {
  r <- assays$rna[[]]; a <- assays$atac[[]]
  cells <- intersect(rownames(r), rownames(a))
  data.frame(sample_id = sample, barcode = cells,
    r[cells, qc_metric_names[1:4], drop = FALSE],
    a[cells, qc_metric_names[5:8], drop = FALSE], check.names = FALSE)
}
read_analysis_config <- function(file, sample = NULL) {
  cfg <- jsonlite::fromJSON(file, simplifyVector = FALSE)
  if (is.null(cfg$thresholds)) stop_bad("QC JSON must have thresholds")
  if (!setequal(names(cfg$thresholds), qc_metric_names))
    stop_bad("QC JSON thresholds must name exactly: ", paste(qc_metric_names, collapse = ", "))
  for (key in c("sample_overrides")) {
    value <- cfg[[key]]
    if (!is.null(value) && (!is.list(value) ||
        (length(value) && (is.null(names(value)) || any(!nzchar(names(value)))))))
      stop_bad(key, " must be an object keyed by sample ID")
    for (entry in value)
      if (!is.list(entry) || is.null(names(entry)) ||
          !all(names(entry) %in% qc_metric_names)) stop_bad("Invalid sample override")
  }
  overrides <- if (!is.null(sample)) cfg$sample_overrides[[sample]] else NULL
  if (!is.null(overrides) && !all(names(overrides) %in% qc_metric_names))
    stop_bad("Unknown override metric for ", sample)
  limits <- cfg$thresholds
  for (nm in names(overrides)) limits[[nm]] <- overrides[[nm]]
  for (nm in qc_metric_names) {
    bounds <- limits[[nm]]
    if (!is.list(bounds) || !setequal(names(bounds), c("min", "max")))
      stop_bad("Expected min and max for ", nm)
    if (all(vapply(bounds, is.null, logical(1)))) stop_bad("No bounds set for ", nm)
    for (v in bounds) if (!is.null(v) && (length(v) != 1L || !is.numeric(v) ||
                                       !is.finite(v))) stop_bad("Invalid bound for ", nm)
    if (!is.null(bounds$min) && !is.null(bounds$max) && bounds$min > bounds$max)
      stop_bad("min exceeds max for ", nm)
  }
  clustering <- cfg$clustering
  # Older configurations did not have an RNA integration switch.
  if (is.list(clustering) && is.null(clustering$integrate_rna)) clustering$integrate_rna <- TRUE
  if (!is.list(clustering) || !setequal(names(clustering),
      c("algorithm", "resolution", "RNA_dims", "RNA_umap_dims", "ATAC_dims",
        "ATAC_integrate_dims", "integrate_rna", "integrate_atac", "min_cutoff")))
    stop_bad("QC JSON clustering must contain algorithm, resolution, RNA_dims, ",
             "RNA_umap_dims, ATAC_dims, ATAC_integrate_dims, integrate_rna, integrate_atac, min_cutoff")
  algorithm <- clustering$algorithm
  resolution <- clustering$resolution
  if (length(algorithm) != 1L || !is.numeric(algorithm) ||
      !is.finite(algorithm) || !(algorithm %in% 1:4))
    stop_bad("clustering.algorithm must be 1, 2, 3, or 4")
  if (length(resolution) != 1L || !is.numeric(resolution) ||
      !is.finite(resolution) || resolution <= 0)
    stop_bad("clustering.resolution must be a positive number")
  if (length(clustering$integrate_atac) != 1L || !is.logical(clustering$integrate_atac) ||
      is.na(clustering$integrate_atac)) stop_bad("integrate_atac must be true or false")
  if (length(clustering$integrate_rna) != 1L || !is.logical(clustering$integrate_rna) ||
      is.na(clustering$integrate_rna)) stop_bad("integrate_rna must be true or false")
  cutoff <- clustering$min_cutoff
  if (length(cutoff) != 1L || !is.numeric(cutoff) ||
      !is.finite(cutoff) || cutoff < 0 || cutoff != floor(cutoff))
    stop_bad("min_cutoff must be a nonnegative integer")
  parse_dims <- function(x, name) {
    if (length(x) != 1L || !is.character(x) ||
        !grepl("^[1-9][0-9]*:[1-9][0-9]*$", x))
      stop_bad(name, " must be an inclusive range such as '1:20'")
    endpoints <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    if (anyNA(endpoints) || endpoints[1] > endpoints[2] ||
        endpoints[2] > 500L) stop_bad("Invalid dimension range for ", name)
    seq.int(endpoints[1], endpoints[2])
  }
  cfg$thresholds <- limits
  cfg$clustering <- list(algorithm = as.integer(algorithm), resolution = resolution,
    RNA_dims = parse_dims(clustering$RNA_dims, "RNA_dims"),
    RNA_umap_dims = parse_dims(clustering$RNA_umap_dims, "RNA_umap_dims"),
    ATAC_dims = parse_dims(clustering$ATAC_dims, "ATAC_dims"),
    ATAC_integrate_dims = parse_dims(clustering$ATAC_integrate_dims, "ATAC_integrate_dims"),
    integrate_rna = clustering$integrate_rna,
    integrate_atac = clustering$integrate_atac, min_cutoff = as.integer(cutoff))
  if (clustering$integrate_atac &&
      max(cfg$clustering$ATAC_dims) > max(cfg$clustering$ATAC_integrate_dims))
    stop_bad("ATAC_integrate_dims must include at least the largest ATAC_dims component")
  doublet <- cfg$doublet
  if (is.null(doublet)) doublet <- list()
  if (!is.list(doublet) || !all(names(doublet) %in%
      c("expected_rate", "sample_rates", "exclude_clusters"))) stop_bad("Invalid doublet settings")
  for (key in c("sample_rates", "exclude_clusters"))
    if (!is.null(doublet[[key]]) && !is.list(doublet[[key]]))
      stop_bad("doublet.", key, " must be an object keyed by sample ID")
  rate <- if (is.null(doublet$expected_rate)) 0.05 else doublet$expected_rate
  if (!is.null(sample) && !is.null(doublet$sample_rates[[sample]]))
    rate <- doublet$sample_rates[[sample]]
  for (v in c(list(rate), doublet$sample_rates))
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v <= 0 || v >= 1)
      stop_bad("Doublet rates must be numbers between zero and one")
  for (v in doublet$exclude_clusters) {
    ids <- unlist(v, use.names = FALSE)
    if (length(ids) && (!is.character(ids) || anyNA(ids) || any(!nzchar(ids))))
      stop_bad("exclude_clusters must contain arrays of cluster ID strings")
  }
  excluded <- if (is.null(sample)) character() else
    as.character(unlist(doublet$exclude_clusters[[sample]], use.names = FALSE))
  cfg$doublet <- list(expected_rate = rate, exclude_clusters = excluded)
  cfg
}
preprocess_qc_rna <- function(obj, settings) {
  set.seed(123)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  pcs <- max(settings$RNA_dims)
  if (min(ncol(obj), length(VariableFeatures(obj))) <= pcs)
    stop_bad("Insufficient cells/features for requested per-sample RNA PCA")
  obj <- RunPCA(obj, npcs = pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = settings$RNA_dims, verbose = FALSE)
  FindClusters(obj, resolution = settings$resolution, algorithm = settings$algorithm,
    cluster.name = "pre_doublet_clusters", verbose = FALSE)
}
remove_qc_clusters <- function(obj, excluded) {
  unknown <- setdiff(excluded, as.character(unique(obj$pre_doublet_clusters)))
  if (length(unknown)) stop_bad("Unknown pre-doublet cluster IDs: ", paste(unknown, collapse = ", "))
  cells <- colnames(obj)[!as.character(obj$pre_doublet_clusters) %in% excluded]
  if (length(cells) < 100L) stop_bad("Fewer than 100 cells after cluster exclusion")
  subset(obj, cells = cells)
}
score_sample_cell_cycle <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  # Match the Seurat 2019 marker panel to the gene symbols present in mouse RNA.
  # This is symbol matching, not a comprehensive human-to-mouse orthology map.
  match_mouse <- function(genes) {
    symbols <- rownames(obj)
    unique(stats::na.omit(symbols[match(toupper(genes), toupper(symbols))]))
  }
  s <- match_mouse(Seurat::cc.genes.updated.2019$s.genes)
  g2m <- match_mouse(Seurat::cc.genes.updated.2019$g2m.genes)
  if (length(s) < 5L || length(g2m) < 5L)
    stop_bad("Too few cell-cycle markers; RNA inputs must use mouse gene symbols")
  # Keep control sampling within the expression-bin size for small gene panels.
  ctrl <- max(1L, min(length(s), length(g2m), floor(nrow(obj) / 24L) - 1L))
  obj <- CellCycleScoring(obj, s.features = s, g2m.features = g2m,
    ctrl = ctrl, set.ident = FALSE)
  obj@misc$cell_cycle_features <- list(S = s, G2M = g2m)
  obj
}
annotate_sample_rna <- function(obj) {
  if (!requireNamespace("cellXY", quietly = TRUE))
    stop_bad("Install phipsonlab/cellXY for classifySex")
  obj <- score_sample_cell_cycle(obj)
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  prediction <- cellXY::classifySex(counts, genome = "Mm", qc = FALSE)
  if (!"prediction" %in% colnames(prediction) ||
      !all(colnames(obj) %in% rownames(prediction)))
    stop_bad("classifySex returned missing cell IDs or prediction column")
  sex <- as.character(prediction[colnames(obj), "prediction"])
  sex[sex == "NA"] <- NA_character_
  obj$sex_prediction <- sex
  obj
}
read_stage_objects <- function(files, assay) {
  items <- lapply(files, readRDS)
  ids <- vapply(items, function(x) {
    ids <- unique(as.character(x$sample_id))
    if (length(ids) != 1L || is.na(ids) || !nzchar(ids))
      stop_bad("Expected one sample_id per input object")
    if (!(assay %in% Assays(x))) stop_bad("Missing ", assay, " assay for ", ids)
    ids
  }, character(1))
  if (anyDuplicated(ids)) stop_bad("Duplicate sample ID")
  if (anyDuplicated(unlist(lapply(items, colnames)))) stop_bad("Duplicate cell ID")
  names(items) <- ids
  items[order(ids)]
}
cluster_qc_reduction <- function(obj, reduction, dims, umap_dims, prefix, settings) {
  if (ncol(Embeddings(obj, reduction)) < max(c(dims, umap_dims)))
    stop_bad("Insufficient dimensions in ", reduction)
  graph <- paste0(prefix, "_snn")
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims,
    graph.name = c(paste0(prefix, "_nn"), graph), verbose = FALSE)
  obj <- FindClusters(obj, graph.name = graph, algorithm = settings$algorithm,
    resolution = settings$resolution, cluster.name = paste0(prefix, "_clusters"),
    verbose = FALSE)
  RunUMAP(obj, reduction = reduction, dims = umap_dims,
    reduction.name = paste0(prefix, "_umap"),
    reduction.key = paste0(gsub("_", "", prefix), "UMAP_"), verbose = FALSE)
}
apply_qc <- function(tab, limits) {
  keep <- rep(TRUE, nrow(tab))
  for (nm in qc_metric_names) {
    v <- tab[[nm]]; b <- limits[[nm]]
    ok <- is.finite(v)
    if (!is.null(b$min)) ok <- ok & v >= b$min
    if (!is.null(b$max)) ok <- ok & v <= b$max
    keep <- keep & ok
  }
  tab$barcode[keep]
}
