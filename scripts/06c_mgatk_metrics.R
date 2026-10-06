#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

  stopifnot(length(args) == 4L)
  sample <- args[1]; obj <- readRDS(args[2]); dir <- args[3]; outfile <- args[4]
  m <- ReadMGATK(dir)
  if (!"A" %in% m$refallele$ref[m$refallele$pos == 5019])
    stop_bad("mgatk reference allele at m.5019 is not A: ", sample)
  bcs <- obj$barcode
  if (!any(bcs %in% colnames(m$counts)))
    stop_bad("No shared barcodes between mgatk and Seurat: ", sample)
  base_counts <- function(base, strand) {
    key <- paste(base, 5019, strand, sep = "-")
    if (!key %in% rownames(m$counts)) stop_bad("Missing mgatk row ", key)
    ans <- rep(NA_real_, length(bcs))
    ix <- match(bcs, colnames(m$counts))
    ok <- which(!is.na(ix))
    ans[ok] <- as.numeric(m$counts[key, ix[ok], drop = TRUE])
    ans
  }
  A_f <- base_counts("A", "fwd"); A_r <- base_counts("A", "rev")
  G_f <- base_counts("G", "fwd"); G_r <- base_counts("G", "rev")
  C_f <- base_counts("C", "fwd"); C_r <- base_counts("C", "rev")
  T_f <- base_counts("T", "fwd"); T_r <- base_counts("T", "rev")
  A <- A_f + A_r; G <- G_f + G_r
  site_depth <- A + G + C_f + C_r + T_f + T_r
  informative <- A + G
  hf <- ifelse(informative > 0, G / informative, NA_real_)
  depth <- m$depth[bcs, 1]
  out <- data.frame(cell_id = colnames(obj), sample_id = sample,
    barcode = bcs, mtDNA_depth = as.numeric(depth),
    depth_5019 = site_depth, informative_depth_5019 = informative,
    A_5019 = A, G_5019 = G, G_5019_fwd = G_f, G_5019_rev = G_r,
    hf_5019_fwd = ifelse(A_f + G_f > 0, G_f / (A_f + G_f), NA_real_),
    hf_5019_rev = ifelse(A_r + G_r > 0, G_r / (A_r + G_r), NA_real_),
    hf_5019 = hf, hf_5019_covered = ifelse(informative >= 10, hf, NA_real_))
  write_tsv(out, outfile)
