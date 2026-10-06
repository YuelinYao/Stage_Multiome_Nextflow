#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
rna <- readRDS(args[1]); atac <- readRDS(args[2])
if (!setequal(colnames(rna), colnames(atac))) stop_bad("RNA/ATAC cells differ")
atac <- atac[, colnames(rna)]
stopifnot(identical(colnames(rna), colnames(atac)))
rna[["ATAC"]] <- atac[["ATAC"]]
rna[["lsi"]] <- atac[["lsi"]]
for (key in c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"))
  rna[[key]] <- atac[[key]][colnames(rna), 1, drop = TRUE]
saveRDS(rna, args[3])
