#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

# Usage: sample_id retained_rna.rds mgatk_enriched_final_dir sample.mtDNA.tsv
# The 5019 and 5024 variant panels follow YY_06_mgatk_heteroplasmy.R.
# Heteroplasmy_<position> and Depth_<position> use its metadata column names.
# This task computes per-cell metrics; step 13 attaches them to the WNN object.
stopifnot(length(args) == 4L)
sample <- args[1]; obj <- readRDS(args[2]); dir <- args[3]; outfile <- args[4]
# ReadMGATK would also match macOS "._" metadata files; read via links to the real files only.
files <- list.files(dir, full.names = TRUE)
files <- files[!startsWith(basename(files), "._")]
clean_dir <- file.path(tempdir(), "mgatk")
dir.create(clean_dir)
if (length(files))
  invisible(file.symlink(normalizePath(files), file.path(clean_dir, basename(files))))
m <- ReadMGATK(clean_dir)
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
