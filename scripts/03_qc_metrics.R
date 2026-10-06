#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "00_common.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
sample <- args[1]; arc <- args[2]; bed <- args[3]
tab <- qc_table(sample, make_assays(arc, bed))
write_tsv(tab, args[4])
plot_qc_violins(tab, args[5], sample)
