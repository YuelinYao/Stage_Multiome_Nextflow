#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

# Usage: sample_id retained_rna.rds mgatk_enriched_dir mgatk_dir sample.mtDNA.tsv [settings.json]
# Either mgatk directory may be "none". Heteroplasmy uses the enriched run when
# supplied, otherwise the ordinary run; coverage is plotted for every supplied run.
# The 5019 and 5024 variant panels follow YY_06_mgatk_heteroplasmy.R.
# Heteroplasmy_<position> and Depth_<position> use its metadata column names.
# This task computes per-cell metrics; step 13 attaches them to the WNN object.
stopifnot(length(args) %in% 5:6)
sample <- args[1]; obj <- readRDS(args[2]); outfile <- args[5]
# Variant thresholds from the MGATK settings section; defaults without a settings file.
variant_settings <- if (length(args) == 6L) read_pipeline_config(args[6], process = "MGATK")$MGATK else
  list(vmr = 0.01, strand_correlation = 0.65, n_cells_conf_detected = 300)
dirs <- c(enriched = args[3], ordinary = args[4])
dirs <- dirs[dirs != "none"]
if (!length(dirs)) stop_bad("Provide an enriched and/or ordinary mgatk directory: ", sample)
# ReadMGATK would also match macOS "._" metadata files; read via links to the real files only.
read_mgatk_clean <- function(dir) {
  files <- list.files(dir, full.names = TRUE)
  files <- files[!startsWith(basename(files), "._")]
  clean_dir <- tempfile("mgatk")
  dir.create(clean_dir)
  if (length(files))
    invisible(file.symlink(normalizePath(files), file.path(clean_dir, basename(files))))
  ReadMGATK(clean_dir)
}
het_source <- names(dirs)[1]
m <- read_mgatk_clean(dirs[[het_source]])
message(sample, ": heteroplasmy from ", het_source, " mgatk ", dirs[[het_source]])
variants <- data.frame(
  pos = c(1781L, 1866L, 3009L, 3823L, 5019L, 5024L, 13614L, 13715L, 15200L, 16232L),
  ref = c("C", "A", "G", "T", "A", "C", "C", "C", "A", "A"),
  alt = c("T", "G", "T", "C", "G", "T", "T", "T", "G", "T"))
# A different reference changes the biological meaning of these allele fractions.
# Check every site, rather than checking only 5019 and silently using the others.
for (i in seq_len(nrow(variants))) {
  ref <- m$refallele$ref[m$refallele$pos == variants$pos[i]]
  if (length(ref) != 1L || is.na(ref) || toupper(ref) != variants$ref[i])
    stop_bad("mgatk reference allele at m.", variants$pos[i], " is not ",
             variants$ref[i], ": ", sample)
}
bcs <- as.character(obj$barcode)
if (length(bcs) != ncol(obj) || anyNA(bcs) || any(!nzchar(bcs)) || anyDuplicated(bcs))
  stop_bad("Missing or duplicated original barcodes: ", sample)
if (is.null(colnames(obj)) || anyNA(colnames(obj)) || anyDuplicated(colnames(obj)))
  stop_bad("Missing or duplicated Seurat cell IDs: ", sample)
for (ids in list(colnames(m$counts), rownames(m$counts), rownames(m$depth)))
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop_bad("Missing or duplicated mgatk identifiers: ", sample)
# Check stored sparse values without expanding the genome-wide count matrix.
values <- if (inherits(m$counts, "sparseMatrix")) m$counts@x else as.numeric(m$counts)
if (any(!is.finite(values) | values < 0 | values != floor(values)))
  stop_bad("mgatk counts must be finite nonnegative integers: ", sample)
ix <- match(bcs, colnames(m$counts))
present <- which(!is.na(ix))
if (!length(present)) stop_bad("No shared barcodes between mgatk and Seurat: ", sample)
message(sample, ": mgatk matches ", length(present), "/", length(bcs), " retained cells")

# Exact row names remove any dependence on mgatk's base/strand row order.
# A zero count for a matched barcode is observed zero coverage; an unmatched
# barcode is missing data. Neither should become a confident wild-type call.
base_counts <- function(pos, base, strand) {
  key <- paste(base, pos, strand, sep = "-")
  if (!key %in% rownames(m$counts)) stop_bad("Missing mgatk row ", key)
  ans <- rep(NA_real_, length(bcs))
  ans[present] <- as.numeric(m$counts[key, ix[present], drop = TRUE])
  ans
}
allele_fraction <- function(alt, ref, min_depth = 1) {
  informative <- alt + ref
  ans <- rep(NA_real_, length(informative))
  usable <- which(is.finite(informative) & informative > 0 & informative >= min_depth)
  ans[usable] <- alt[usable] / informative[usable]
  ans
}
if (ncol(m$depth) != 1L) stop_bad("Expected one mgatk depth column: ", sample)
depth <- as.numeric(m$depth[match(bcs, rownames(m$depth)), 1])
if (any(!is.na(depth) & (!is.finite(depth) | depth < 0)))
  stop_bad("Invalid mgatk sequencing depth: ", sample)
depth[is.na(ix)] <- NA_real_
mt_count <- rep(NA_real_, length(bcs))
mt_count[present] <- as.numeric(Matrix::colSums(m$counts[, ix[present], drop = FALSE]))
out <- data.frame(cell_id = colnames(obj), sample_id = sample, barcode = bcs,
                  mgatk_matched = !is.na(ix), mtDNA_depth = depth, mtDNA_count = mt_count)
sites <- list()
for (i in seq_len(nrow(variants))) {
  pos <- variants$pos[i]
  counts <- lapply(c("A", "C", "G", "T"), function(base)
    list(fwd = base_counts(pos, base, "fwd"), rev = base_counts(pos, base, "rev")))
  names(counts) <- c("A", "C", "G", "T")
  ref <- counts[[variants$ref[i]]]; alt <- counts[[variants$alt[i]]]
  ref_total <- ref$fwd + ref$rev; alt_total <- alt$fwd + alt$rev
  site <- list(depth = Reduce(`+`, lapply(counts, function(x) x$fwd + x$rev)),
               ref = ref_total, alt = alt_total)
  sites[[as.character(pos)]] <- site
  out[[paste0("Depth_", pos)]] <- site$depth
  out[[paste0("Informative_depth_", pos)]] <- ref_total + alt_total
  out[[paste0("Allele_ref_", pos)]] <- ref_total
  out[[paste0("Allele_mut_", pos)]] <- alt_total
  # Export all eight counts, including other bases, for auditing strand balance.
  for (base in names(counts)) {
    out[[paste0(base, "_", pos)]] <- counts[[base]]$fwd + counts[[base]]$rev
    for (strand in c("fwd", "rev"))
      out[[paste(base, pos, strand, sep = "_")]] <- counts[[base]][[strand]]
  }
  out[[paste0("Heteroplasmy_", pos)]] <- allele_fraction(alt_total, ref_total)
  out[[paste0("Heteroplasmy_", pos, "_covered")]] <- allele_fraction(alt_total, ref_total, 10)
  out[[paste0("Heteroplasmy_", pos, "_fwd")]] <- allele_fraction(alt$fwd, ref$fwd)
  out[[paste0("Heteroplasmy_", pos, "_rev")]] <- allele_fraction(alt$rev, ref$rev)
}

# Pool READ COUNTS, not site fractions: sites with more informative reads have
# more weight. This assumes the supplied sites tag the same haplotype; these
# short-read summaries do not establish linkage or phase. Other-base reads
# contribute to total depth but not to the two-allele fraction or its cutoff.
add_pool <- function(label, positions, reciprocal = integer()) {
  # Reciprocal sites exchange reference/alternate roles only in this pool.
  # Individual-site columns always retain the reference-genome orientation.
  stopifnot(all(reciprocal %in% positions))
  selected <- sites[as.character(positions)]
  for (pos in reciprocal) {
    key <- as.character(pos)
    original <- selected[[key]]
    selected[[key]]$ref <- original$alt
    selected[[key]]$alt <- original$ref
  }
  ref <- Reduce(`+`, lapply(selected, `[[`, "ref"))
  alt <- Reduce(`+`, lapply(selected, `[[`, "alt"))
  out[[paste0("Depth_", label)]] <<- Reduce(`+`, lapply(selected, `[[`, "depth"))
  out[[paste0("Informative_depth_", label)]] <<- ref + alt
  out[[paste0("Allele_ref_", label)]] <<- ref
  out[[paste0("Allele_mut_", label)]] <<- alt
  out[[paste0("Heteroplasmy_", label)]] <<- allele_fraction(alt, ref)
  # YY_06 displays pooled results at depth >=25. Require >=25 informative
  # reads here, so unrelated base calls cannot satisfy the confidence cutoff.
  out[[paste0("Heteroplasmy_", label, "_covered")]] <<- allele_fraction(alt, ref, 25)
}
add_pool("all_5019", c(5019L, 15200L, 16232L))
# Apply the same 5024 haplotype model to every sample (including 5019 and WT).
# 5024 and 13715 alternative alleles tag the 5024-mutant mitochondria.
# At the other five sites, alternative alleles tag mitochondria with WT 5024;
# reverse those sites so their reference counts contribute to the mutant pool.
# Do not branch on sample names or force WT/5019 measurements to zero.
positions_5024 <- c(1781L, 3009L, 5024L, 13614L, 13715L, 3823L, 1866L)
add_pool("all_5024", positions_5024,
         reciprocal = c(1781L, 3009L, 13614L, 3823L, 1866L))
write_tsv(out, outfile)

# chrM coverage: mean per retained cell at each position, for every supplied
# mgatk run. A+C+G+T counts on both strands equal mgatk's own coverage values.
mean_coverage <- function(mg) {
  cols <- match(bcs, colnames(mg$counts))
  cols <- cols[!is.na(cols)]
  pos <- as.integer(sub("^[ACGT]-([0-9]+)-.*$", "\\1", rownames(mg$counts)))
  per_row <- if (length(cols)) Matrix::rowSums(mg$counts[, cols, drop = FALSE]) / length(cols) else 0 * pos
  sums <- rowsum(per_row, pos)
  setNames(sums[, 1], rownames(sums))
}
coverage <- data.frame(position = sort(unique(m$refallele$pos)))
for (src in names(dirs)) {
  mg <- if (identical(dirs[[src]], dirs[[het_source]])) m else read_mgatk_clean(dirs[[src]])
  values <- unname(mean_coverage(mg)[as.character(coverage$position)])
  coverage[[paste0("mean_coverage_", src)]] <- ifelse(is.na(values), 0, values)
}
prefix <- file.path(dirname(outfile), paste0(sample, ".mtDNA_coverage"))
write_tsv(coverage, paste0(prefix, ".tsv"))
# Enriched (green) is drawn first and ordinary (yellow) on top; means <= 1 are left blank.
long <- do.call(rbind, lapply(names(dirs), function(src) data.frame(position = coverage$position,
  source = src, mean_coverage = coverage[[paste0("mean_coverage_", src)]])))
long <- long[long$mean_coverage > 1, , drop = FALSE]
long$source <- factor(long$source, levels = c("enriched", "ordinary"))
p <- ggplot2::ggplot(long, ggplot2::aes(x = position, y = mean_coverage, fill = source)) +
  ggplot2::geom_col(width = 1, position = "identity") +
  ggplot2::scale_fill_manual(values = c(enriched = "#64b53b", ordinary = "#fdcb25")) +
  ggplot2::scale_y_continuous(trans = "log10") +
  ggplot2::labs(x = "Position along chrM", y = "Mean coverage per cell (log10)",
                title = sample, fill = "mgatk") +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(1, 1, 1, 1))
ggplot2::ggsave(paste0(prefix, ".pdf"), plot = p, width = 4, height = 3)
# Saved ggplot object, for combining samples outside the pipeline.
saveRDS(p, paste0(prefix, "_plot.rds"))

# Variant discovery with Signac IdentifyVariants on the retained cells of the
# heteroplasmy mgatk run. High-confidence variants pass all three MGATK thresholds;
# the plot shows variants with at least n_cells_conf_detected confident cells.
cells <- colnames(m$counts)[ix[present]]
variants <- identify_mtdna_variants(m$counts[, cells, drop = FALSE], m$refallele)
variants$high_confidence <- with(variants,
  !is.na(vmr) & !is.na(strand_correlation) &
  n_cells_conf_detected >= variant_settings$n_cells_conf_detected &
  strand_correlation > variant_settings$strand_correlation &
  vmr > variant_settings$vmr)
message(sample, ": IdentifyVariants found ", nrow(variants), " variants, ",
        sum(variants$high_confidence), " high-confidence")
variant_prefix <- file.path(dirname(outfile), paste0(sample, ".mtDNA_variants"))
write_tsv(variants, paste0(variant_prefix, ".tsv"))
write_tsv(variants[variants$high_confidence, , drop = FALSE], paste0(variant_prefix, "_high_conf.tsv"))
plot_df <- variants[variants$n_cells_conf_detected >= variant_settings$n_cells_conf_detected &
                    is.finite(log10(variants$vmr)) & !is.na(variants$strand_correlation), , drop = FALSE]
p_variants <- ggplot2::ggplot(plot_df, ggplot2::aes(x = strand_correlation, y = log10(vmr),
                                                    color = high_confidence)) +
  ggplot2::geom_point(size = 0.4) +
  ggplot2::scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "firebrick")) +
  ggplot2::geom_vline(xintercept = variant_settings$strand_correlation, linetype = 2) +
  ggplot2::geom_hline(yintercept = log10(variant_settings$vmr), linetype = 2) +
  ggrepel::geom_text_repel(data = plot_df[plot_df$high_confidence, , drop = FALSE],
                           ggplot2::aes(label = variant), size = 2.3, max.overlaps = Inf) +
  ggplot2::scale_x_continuous(breaks = c(0, 0.7, 0.8, 0.9)) +
  ggplot2::labs(x = "Strand concordance", y = "log VMR", title = sample) +
  BuenColors::pretty_plot(fontsize = 7.5) +
  BuenColors::L_border() +
  ggplot2::theme(legend.position = "none", axis.text = ggplot2::element_text(size = 9))
ggplot2::ggsave(paste0(variant_prefix, ".pdf"), plot = p_variants, width = 3, height = 3)
saveRDS(p_variants, paste0(variant_prefix, "_plot.rds"))
