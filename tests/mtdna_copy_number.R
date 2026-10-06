# Synthetic integration tests. Run from the repository root with:
# Rscript tests/mtdna_copy_number.R
source("scripts/11_estimate_mtdna_copy_number.R")
for (pkg in c("Rsamtools", "GenomicRanges", "GenomicAlignments"))
  stopifnot(requireNamespace(pkg, quietly = TRUE))

# Test CIGAR-aware base counting: insertions/deletions/soft clips do not create
# nuclear reference coverage; low-quality and ambiguous bases are excluded.
intervals <- IRanges::IRanges(c(1, 8), c(3, 10))
stopifnot(count_read_bases("5M2D5M", 1L, "AAAAAAAAAA", "IIIIIIIIII", intervals, 20) == 6)
stopifnot(count_read_bases("5M2I5M", 1L, "AAAAAAAAAAAA", "IIIIIIIIIIII", intervals, 20) == 6)
stopifnot(count_read_bases("2S10M", 1L, "AAAAAAAAAAAA", "IIIIIIIIIIII", intervals, 20) == 6)
stopifnot(count_read_bases("3M", 1L, "ANA", "I!I", intervals, 20) == 2)
stopifnot(count_read_bases("3M", 1L, "AAA", "!!!", intervals, 20) == 0)
stopifnot(count_read_bases("3M", 1L, "AAA", NA_character_, intervals, 20) == 0,
          count_read_bases("3M", 1L, "*", "*", intervals, 20) == 0)

# A tiny coordinate-sorted BAM exercises real BAM reading, CB translation,
# exclusion of duplicates/secondary alignments, quality filtering, and a read
# spanning two disjoint reference intervals (counted once, only inside BED).
run_tests <- function() {
  root <- tempfile("mtdna-test-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  sam <- file.path(root, "nuclear.sam")
  alignment <- function(name, flag, pos, cigar, sequence, quality, cb, mapq = 60)
    paste(name, flag, "chr1", pos, mapq, cigar, "*", 0, 0, sequence, quality,
          paste0("CB:Z:", cb), sep = "\t")
  writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:2000",
    alignment("usable", 0, 1, "10M", "AAAAAAAAAA", "IIIIIIIIII", "ATAC1"),
    alignment("duplicate", 1024, 1, "10M", "AAAAAAAAAA", "IIIIIIIIII", "ATAC1"),
    alignment("secondary", 256, 1, "10M", "AAAAAAAAAA", "IIIIIIIIII", "ATAC1"),
    alignment("lowmapq", 0, 1, "10M", "AAAAAAAAAA", "IIIIIIIIII", "ATAC1", 10),
    alignment("lowbaseq", 0, 1, "10M", "AAAAAAAAAA", "!!!!!!!!!!", "ATAC1"),
    alignment("missingmt", 0, 1, "10M", "AAAAAAAAAA", "IIIIIIIIII", "CELL2")), sam)
  bam <- Rsamtools::asBam(sam, destination = file.path(root, "nuclear"), overwrite = TRUE)
  bed <- file.path(root, "reference.bed")
  writeLines(c("chr1\t0\t3", "chr1\t7\t10"), bed)
  cells <- file.path(root, "barcodes.txt"); writeLines(c("CELL1", "CELL2", "CELL3"), cells)
  map <- file.path(root, "per_barcode_metrics.csv")
  writeLines(c("barcode,atac_barcode", "CELL1,ATAC1", "CELL2,ATAC2", "CELL3,ATAC3"), map)
  mgatk <- file.path(root, "mgatk"); dir.create(mgatk)
  writeLines(c("ATAC1\t100", "CELL3\t50"), file.path(mgatk, "sample.depthTable.txt"))
  prefix <- file.path(root, "output", "sample")
  args <- c("--sample", "test", "--bam", bam, "--barcodes", cells,
    "--mgatk-dir", mgatk, "--barcode-map", map, "--reference-bed", bed,
    "--min-nuclear-bases", "1", "--out-prefix", prefix)
  out <- main(args)
  stopifnot(out$nuclear_reference_bp[1] == 6,
    out$nuclear_aligned_bases[1] == 6, out$nuclear_reads_with_counted_bases[1] == 1,
    out$nuclear_mean_depth[1] == 1, out$mtDNA_copy_estimate[1] == 200,
    out$qc_status[2] == "missing_mgatk", is.na(out$mtDNA_copy_estimate[2]),
    out$qc_status[3] == "zero_nuclear_coverage", is.na(out$mtDNA_copy_estimate[3]))
  stopifnot(file.exists(paste0(prefix, ".settings.txt")),
    identical(read_bed(paste0(prefix, ".nuclear_reference.bed"), c(chr1 = 2000)),
              read_bed(bed, c(chr1 = 2000))))
  # Low coverage is flagged separately from missing measurements. Raw ratios
  # are retained, while the qualified estimate is suppressed.
  low_args <- args
  low_args[match("--min-nuclear-bases", low_args) + 1L] <- "1000"
  low <- main(low_args)
  stopifnot(low$qc_status[1] == "low_nuclear_coverage",
            low$mtDNA_copy_estimate_raw[1] == 200,
            is.na(low$mtDNA_copy_estimate[1]))
  # A duplicate-retaining run changes coverage in the expected direction.
  duplicated <- main(c(args, "--keep-duplicates", "true"))
  stopifnot(duplicated$nuclear_aligned_bases[1] == 12,
            duplicated$mtDNA_copy_estimate[1] == 100)
  # Fragment endpoints are scored once each, ignoring duplicate readSupport.
  # A fully blacklisted bin cannot enter selection. chrM is never eligible.
  fragments <- file.path(root, "fragments.tsv.gz")
  con <- gzfile(fragments, "wt")
  writeLines(c("# test", "chr1\t0\t10\tCELL1\t100", "chr1\t0\t11\tCELL1\t1",
               "chr1\t500\t510\tCELL1\t1", "chrM\t0\t10\tCELL1\t100"), con)
  close(con)
  selected <- select_reference(fragments, c("CELL1"), c(chr1 = 2000, chrM = 100),
    GenomicRanges::GRanges("chr1", IRanges::IRanges(501, 1000)), 500, 0.5)
  stopifnot(length(selected) == 1, IRanges::start(selected) == 1,
            IRanges::end(selected) == 500)
  # Full automatic-reference route with ATAC aliases in fragments as well as
  # BAM and mgatk. Before the fix, canonical CELL1 never matched these fragments.
  con <- gzfile(fragments, "wt")
  writeLines(c("chr1\t0\t10\tATAC1\t100", "chr1\t0\t11\tATAC1\t1",
               "chr1\t500\t510\tATAC1\t1"), con)
  close(con)
  blacklist <- file.path(root, "blacklist.bed")
  writeLines("chr1\t500\t1000", blacklist)
  selected_args <- args[-c(match("--reference-bed", args), match("--reference-bed", args) + 1L)]
  selected_args <- c(selected_args, "--fragments", fragments, "--blacklist", blacklist,
                     "--quantile", "0.5")
  automatic <- main(selected_args)
  stopifnot(automatic$nuclear_reference_bp[1] == 500,
            automatic$nuclear_aligned_bases[1] == 10,
            automatic$nuclear_mean_depth[1] == 0.02,
            automatic$mtDNA_copy_estimate[1] == 10000)
  expect_error <- function(expr) {
    stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
  }
  writeLines(c("barcode,atac_barcode", "CELL1,ATAC1", "CELL1,ATAC4"), map)
  expect_error(main(args))
  writeLines(c("barcode,atac_barcode", "CELL1,"), map)
  expect_error(main(args))
}
run_tests()
cat("All mtDNA copy-number tests passed.\n")
