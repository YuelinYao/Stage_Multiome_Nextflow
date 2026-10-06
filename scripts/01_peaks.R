#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)

  stopifnot(length(args) == 4L)
  sample <- args[1]; arc <- args[2]; macs <- args[3]; outfile <- args[4]
  x <- read_arc(arc)
  frag <- file.path(arc, "atac_fragments.tsv.gz")
  atac <- CreateChromatinAssay(counts = x$Peaks, sep = c(":", "-"), fragments = frag)
  obj <- CreateSeuratObject(atac, assay = "ATAC", project = sample)
  # Signac 1.x splits the fragment file when group.by is set. Without this,
  # CallPeaks can use all fragments, including non-cell barcodes.
  obj$peak_group <- "filtered_cells"
  peaks <- CallPeaks(obj, assay = "ATAC", group.by = "peak_group", macs2.path = macs,
                     effective.genome.size = 1.87e9, name = sample)
  peaks <- peaks[grepl("^chr([0-9]+|X|Y)$", as.character(GenomicRanges::seqnames(peaks)))]
  if (length(peaks) == 0) stop_bad("No nuclear MACS peaks: ", sample)
  rtracklayer::export(peaks, outfile, format = "BED")
