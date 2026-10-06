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
# Signac 1.14 (the exported environment) does not accept keep_all_features.
# Older FeatureMatrix versions omit peaks on absent chromosomes. Restore those
# rows as sparse zeros, and impose the common BED's exact row/cell order.
count_common_peaks <- function(fragments, peaks, cells) {
  features <- paste(as.character(GenomicRanges::seqnames(peaks)),
    IRanges::start(peaks), IRanges::end(peaks), sep = "-")
  if (!length(features) || anyDuplicated(features) || !length(cells) || anyDuplicated(cells))
    stop_bad("Common peaks and cells must be nonempty and unique")
  arguments <- list(fragments = fragments, features = peaks, cells = cells, sep = c("-", "-"))
  if ("keep_all_features" %in% names(formals(FeatureMatrix))) arguments$keep_all_features <- TRUE
  mat <- do.call(FeatureMatrix, arguments)
  if (is.null(rownames(mat)) || is.null(colnames(mat)) ||
      anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat)) ||
      any(!rownames(mat) %in% features) || any(!colnames(mat) %in% cells))
    stop_bad("Unexpected peak/cell identifiers returned by FeatureMatrix")
  # Sparse coordinate remapping preserves counts without expanding the matrix.
  entries <- Matrix::summary(methods::as(mat, "dgCMatrix"))
  Matrix::sparseMatrix(i = match(rownames(mat), features)[entries$i],
    j = match(colnames(mat), cells)[entries$j], x = entries$x,
    dims = c(length(features), length(cells)), dimnames = list(features, cells))
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
  mat <- count_common_peaks(fragments, peaks, cells)
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
# Independent defaults for each analysis step, used by the configuration reader.
pipeline_analysis_defaults <- function() list(
  review_QC_cluster = list(RNA_dims = "1:20", algorithm = 3L, resolution = 2),
  doublet = list(PC_dims = "1:20", algorithm = 3L, resolution = 2),
  RNA_processing = list(RNA_dims = "1:20", RNA_umap_dims = "1:50",
                        algorithm = 3L, resolution = 2, integrate_rna = TRUE),
  ATAC_processing = list(ATAC_dims = "2:40", ATAC_integrate_dims = "1:40",
                         algorithm = 3L, resolution = 2, integrate_atac = TRUE,
                         min_cutoff = "q5"),
  WNN = list(RNA_dims = "1:20", ATAC_dims = "2:40", algorithm = 3L, resolution = 2))
read_pipeline_config <- function(file, sample = NULL) {
  cfg <- jsonlite::fromJSON(file, simplifyVector = FALSE)
  named_object <- function(value, name) {
    if (!is.list(value) || (length(value) &&
        (is.null(names(value)) || anyNA(names(value)) ||
         any(!nzchar(names(value))) || anyDuplicated(names(value)))))
      stop_bad(name, " must be an object with unique, nonempty keys")
  }
  named_object(cfg, "Pipeline JSON")
  allowed <- c("peaks", "note", "thresholds", "sample_overrides", "review_QC_cluster",
               "doublet", "RNA_processing", "ATAC_processing", "WNN", "clustering")
  unknown <- setdiff(names(cfg), allowed)
  if (length(unknown)) stop_bad("Unknown pipeline JSON sections: ", paste(unknown, collapse = ", "))
  peaks <- cfg$peaks
  if (is.null(peaks)) peaks <- list()
  named_object(peaks, "peaks")
  if (!is.list(peaks) || !all(names(peaks) %in% c("peak_min", "peak_max")))
    stop_bad("Pipeline JSON peaks must contain peak_min and/or peak_max")
  if (is.null(peaks$peak_min)) peaks$peak_min <- 20
  if (is.null(peaks$peak_max)) peaks$peak_max <- 10000
  for (value in peaks)
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value < 0 || value != floor(value))
      stop_bad("Peak widths must be nonnegative integers")
  if (peaks$peak_min >= peaks$peak_max) stop_bad("peak_min must be less than peak_max")
  cfg$peaks <- peaks
  if (is.null(cfg$thresholds)) stop_bad("Pipeline JSON must have thresholds")
  named_object(cfg$thresholds, "thresholds")
  if (!setequal(names(cfg$thresholds), qc_metric_names))
    stop_bad("Pipeline JSON thresholds must name exactly: ", paste(qc_metric_names, collapse = ", "))
  for (key in c("sample_overrides")) {
    value <- cfg[[key]]
    if (!is.null(value)) named_object(value, key)
    if (!is.null(value) && (!is.list(value) ||
        (length(value) && (is.null(names(value)) || any(!nzchar(names(value)))))))
      stop_bad(key, " must be an object keyed by sample ID")
    for (entry in value) {
      named_object(entry, "sample override")
      if (!is.list(entry) || is.null(names(entry)) ||
          !all(names(entry) %in% qc_metric_names)) stop_bad("Invalid sample override")
    }
  }
  overrides <- if (!is.null(sample)) cfg$sample_overrides[[sample]] else NULL
  if (!is.null(overrides) && !all(names(overrides) %in% qc_metric_names))
    stop_bad("Unknown override metric for ", sample)
  limits <- cfg$thresholds
  for (nm in names(overrides)) limits[[nm]] <- overrides[[nm]]
  validate_bounds <- function(bounds, nm) {
    named_object(bounds, paste0("Bounds for ", nm))
    if (!is.list(bounds) || !setequal(names(bounds), c("min", "max")))
      stop_bad("Expected min and max for ", nm)
    if (all(vapply(bounds, is.null, logical(1)))) stop_bad("No bounds set for ", nm)
    for (v in bounds) if (!is.null(v) && (length(v) != 1L || !is.numeric(v) ||
                                       !is.finite(v))) stop_bad("Invalid bound for ", nm)
    if (!is.null(bounds$min) && !is.null(bounds$max) && bounds$min > bounds$max)
      stop_bad("min exceeds max for ", nm)
  }
  # Validate shared and all sample-specific settings, including overrides for
  # samples other than the one currently processed, before any filtering.
  for (nm in qc_metric_names) validate_bounds(cfg$thresholds[[nm]], nm)
  for (id in names(cfg$sample_overrides))
    for (nm in names(cfg$sample_overrides[[id]]))
      validate_bounds(cfg$sample_overrides[[id]][[nm]], paste(id, nm, sep = "/"))
  if (!is.null(cfg$clustering))
    stop_bad("The shared clustering section has been replaced by review_QC_cluster, ",
             "doublet, RNA_processing, ATAC_processing and WNN; see examples/pipeline_settings.example.json")
  parse_dims <- function(x, name) {
    if (length(x) != 1L || !is.character(x) ||
        !grepl("^[1-9][0-9]*:[1-9][0-9]*$", x))
      stop_bad(name, " must be an inclusive range such as '1:20'")
    endpoints <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    if (anyNA(endpoints) || endpoints[1] > endpoints[2] || endpoints[2] > 500L)
      stop_bad("Invalid dimension range for ", name)
    seq.int(endpoints[1], endpoints[2])
  }
  defaults <- pipeline_analysis_defaults()
  parse_section <- function(name, extra = character()) {
    supplied <- cfg[[name]]
    if (is.null(supplied)) supplied <- list()
    if (!is.list(supplied) || (length(supplied) &&
        (is.null(names(supplied)) || anyDuplicated(names(supplied)))) ||
        !all(names(supplied) %in% c(names(defaults[[name]]), extra)))
      stop_bad("Invalid settings in ", name)
    result <- defaults[[name]]
    for (key in names(supplied)) result[key] <- supplied[key]
    if (length(result$algorithm) != 1L || !is.numeric(result$algorithm) ||
        !is.finite(result$algorithm) || !(result$algorithm %in% 1:4))
      stop_bad(name, ".algorithm must be 1, 2, 3, or 4")
    result$algorithm <- as.integer(result$algorithm)
    if (length(result$resolution) != 1L || !is.numeric(result$resolution) ||
        !is.finite(result$resolution) || result$resolution <= 0)
      stop_bad(name, ".resolution must be a positive number")
    for (key in intersect(names(result), c("RNA_dims", "RNA_umap_dims", "ATAC_dims",
                                         "ATAC_integrate_dims")))
      result[[key]] <- parse_dims(result[[key]], paste0(name, ".", key))
    for (key in intersect(names(result), c("integrate_rna", "integrate_atac")))
      if (length(result[[key]]) != 1L || !is.logical(result[[key]]) || is.na(result[[key]]))
        stop_bad(name, ".", key, " must be true or false")
    if ("min_cutoff" %in% names(result)) {
      cutoff <- result$min_cutoff
      numeric_cutoff <- length(cutoff) == 1L && is.numeric(cutoff) &&
        is.finite(cutoff) && cutoff >= 0 && cutoff == floor(cutoff)
      percentile_cutoff <- length(cutoff) == 1L && is.character(cutoff) &&
        !is.na(cutoff) && grepl("^q[0-9]+([.][0-9]+)?$", cutoff) &&
        as.numeric(sub("^q", "", cutoff)) < 100
      if (!numeric_cutoff && !percentile_cutoff)
        stop_bad(name, ".min_cutoff must be a nonnegative integer or a percentile ",
                 "string from q0 to below q100, such as 'q5'")
      # Keep percentile strings intact for Signac::FindTopFeatures().
      result$min_cutoff <- cutoff
    }
    result
  }
  cfg$thresholds <- limits
  for (name in c("review_QC_cluster", "RNA_processing", "ATAC_processing", "WNN"))
    cfg[[name]] <- parse_section(name)
  # IntegrateEmbeddings renumbers output columns sequentially; including LSI1
  # preserves the meaning of the original LSI indices used downstream.
  if (cfg$ATAC_processing$integrate_atac &&
      min(cfg$ATAC_processing$ATAC_integrate_dims) != 1L)
    stop_bad("ATAC_processing.ATAC_integrate_dims must start at 1 when integration is enabled")
  if (cfg$ATAC_processing$integrate_atac &&
      !all(c(cfg$ATAC_processing$ATAC_dims, cfg$WNN$ATAC_dims) %in%
           cfg$ATAC_processing$ATAC_integrate_dims))
    stop_bad("ATAC_processing.ATAC_integrate_dims must cover ATAC_processing.ATAC_dims ",
             "and WNN.ATAC_dims when ATAC integration is enabled")
  doublet <- parse_section("doublet", c("expected_rate", "sample_rates", "exclude_clusters"))
  doublet_dims <- parse_dims(if (is.null(doublet$PC_dims)) "1:20" else doublet$PC_dims,
                            "doublet.PC_dims")
  # DoubletFinder recomputes length(PCs) components internally, so require
  # a contiguous range starting at PC1; our default PCA computes 50 PCs.
  if (min(doublet_dims) != 1L || max(doublet_dims) > 50L)
    stop_bad("doublet.PC_dims must start at 1 and end at or below 50")
  for (key in c("sample_rates", "exclude_clusters"))
    if (!is.null(doublet[[key]])) named_object(doublet[[key]], paste0("doublet.", key))
  rate <- if (is.null(doublet$expected_rate)) "auto" else doublet$expected_rate
  if (!is.null(sample) && !is.null(doublet$sample_rates[[sample]]))
    rate <- doublet$sample_rates[[sample]]
  for (v in c(list(if (is.null(doublet$expected_rate)) "auto" else doublet$expected_rate),
              doublet$sample_rates))
    if (!(identical(v, "auto") ||
          (length(v) == 1L && is.numeric(v) && is.finite(v) && v > 0 && v < 1)))
      stop_bad("Doublet rates must be 'auto' or numbers between zero and one")
  for (v in doublet$exclude_clusters) {
    ids <- unlist(v, use.names = FALSE)
    if (length(ids) && (!is.character(ids) || anyNA(ids) || any(!nzchar(ids))))
      stop_bad("exclude_clusters must contain arrays of cluster ID strings")
  }
  excluded <- if (is.null(sample)) character() else
    as.character(unlist(doublet$exclude_clusters[[sample]], use.names = FALSE))
  cfg$doublet <- list(PC_dims = doublet_dims, algorithm = doublet$algorithm,
    resolution = doublet$resolution, expected_rate = rate, exclude_clusters = excluded)
  cfg
}
read_arc_detected_cells <- function(arc) {
  file <- file.path(arc, "summary.csv")
  summary <- read.csv(file, check.names = FALSE, colClasses = "character")
  column <- "Estimated number of cells"
  if (nrow(summary) != 1L || sum(names(summary) == column) != 1L)
    stop_bad("Expected one row and one 'Estimated number of cells' column in ", file)
  count <- suppressWarnings(as.numeric(gsub(",", "", trimws(summary[[column]]))))
  if (length(count) != 1L || !is.finite(count) || count < 1 || count != floor(count))
    stop_bad("Invalid Estimated number of cells in ", file)
  count
}
nextgem_multiome_doublet_rate <- function(detected_cells) {
  # Next GEM Multiome CG000338 Rev F, p22. Interpolate the published
  # recovered-nuclei table; this estimates a multiplet rate, not a measured rate.
  # https://cdn.10xgenomics.com/image/upload/v1666737555/support-documents/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevF.pdf
  # page 22
  recovered <- c(500, seq(1000, 10000, by = 1000))
  rates <- c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046,
             0.054, 0.062, 0.069, 0.077)
  if (length(detected_cells) != 1L || !is.finite(detected_cells) ||
      detected_cells < 500 || detected_cells > 10000)
    stop_bad("Automatic Next GEM Multiome rate requires 500-10000 ARC detected cells; ",
             "set an explicit doublet.expected_rate or sample_rates override")
  approx(recovered, rates, xout = detected_cells)$y
}
preprocess_qc_rna <- function(obj, settings, default_pca = FALSE) {
  set.seed(123)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  pcs <- if (default_pca) 50L else max(settings$RNA_dims)
  if (min(ncol(obj), length(VariableFeatures(obj))) <= pcs)
    stop_bad("Insufficient cells/features for requested per-sample RNA PCA")
  obj <- if (default_pca) RunPCA(obj, verbose = FALSE) else
    RunPCA(obj, npcs = pcs, verbose = FALSE)
  if (ncol(Embeddings(obj, "pca")) < max(settings$RNA_dims))
    stop_bad("Insufficient computed PCs for requested per-sample RNA dimensions")
  obj <- FindNeighbors(obj, dims = settings$RNA_dims, verbose = FALSE)
  FindClusters(obj, resolution = settings$resolution, algorithm = settings$algorithm,
    cluster.name = "pre_doublet_clusters", verbose = FALSE)
}
# Load the exact cluster labels inspected in phase cluster_review. Match by
# barcode, never row position, and reject missing/extra cells or stale metrics.
read_reviewed_qc_clusters <- function(file, sample, tab, keep) {
  review <- read.delim(file, check.names = FALSE,
    colClasses = c(sample_id = "character", barcode = "character",
                   pre_doublet_cluster = "character"))
  required <- c("sample_id", "barcode", "pre_doublet_cluster", qc_metric_names)
  if (!all(required %in% names(review)) || anyDuplicated(names(review)))
    stop_bad("Missing or duplicated columns in cluster review table: ", file)
  if (!nrow(review) || anyNA(review[, required[1:3]]) ||
      any(review$sample_id != sample) || any(!nzchar(review$barcode)) ||
      any(!nzchar(review$pre_doublet_cluster)) || anyDuplicated(review$barcode))
    stop_bad("Invalid sample, barcode or cluster labels in review table: ", file)
  if (!setequal(review$barcode, keep))
    stop_bad("Reviewed cells differ from current QC-passing cells for ", sample,
             "; rerun --phase cluster_review with the current QC thresholds")
  review <- review[match(keep, review$barcode), , drop = FALSE]
  current <- tab[match(keep, tab$barcode), , drop = FALSE]
  for (metric in qc_metric_names)
    if (!is.numeric(review[[metric]]) ||
        !isTRUE(all.equal(review[[metric]], current[[metric]], tolerance = 1e-8,
                         check.attributes = FALSE)))
      stop_bad("Reviewed QC metrics differ for ", sample, "/", metric,
               "; rerun --phase cluster_review with the current inputs")
  setNames(review$pre_doublet_cluster, review$barcode)
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
