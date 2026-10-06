#!/usr/bin/env Rscript
# Estimate mtDNA abundance relative to accessible, diploid nuclear DNA.
#
# Formula: estimated copies = nuclear ploidy * mean mtDNA depth / mean nuclear depth.
# This is an assay-dependent estimate, not an independently calibrated molecule count.
# The script is standalone: it is NOT automatically run by main.nf.
# Run with --help for arguments. All BED files use 0-based, end-exclusive coordinates.
# Required R packages are already listed in environment.yml; Seurat/Signac are not needed.

fail <- function(...) stop(..., call. = FALSE)

# Import BED explicitly: converting start+1 is essential because GRanges uses
# 1-based inclusive coordinates, whereas BED uses 0-based exclusive ends.
read_bed <- function(path, lengths) {
  tab <- read.table(path, sep = "\t", header = FALSE, comment.char = "#",
                    stringsAsFactors = FALSE, quote = "")
  if (ncol(tab) < 3L || !nrow(tab)) fail("Empty/invalid BED: ", path)
  tab <- tab[, 1:3]; names(tab) <- c("chrom", "start", "end")
  if (!is.numeric(tab$start) || !is.numeric(tab$end) ||
      anyNA(tab) || any(tab$start != floor(tab$start)) ||
      any(tab$end != floor(tab$end)) || any(tab$start < 0 | tab$end <= tab$start) ||
      any(!tab$chrom %in% names(lengths)) || any(tab$end > lengths[tab$chrom]))
    fail("BED coordinates/contigs disagree with BAM header: ", path)
  GenomicRanges::GRanges(tab$chrom, IRanges::IRanges(tab$start + 1, tab$end))
}

# Score only full-width autosomal bins. Sex chromosomes are excluded so a
# diploid normalization does not depend on inferred sex. The chromosome lengths
# come from YOUR BAM header, rather than a hard-coded human or mouse genome.
select_reference <- function(fragments, cells, lengths, blacklist, bin_size, quantile_cutoff) {
  chroms <- grep("^chr([1-9]|1[0-9])$", names(lengths), value = TRUE)
  if (!length(chroms)) fail("No mouse chr1-chr19 contigs; provide --reference-bed instead")
  bins <- do.call(c, lapply(chroms, function(chrom) {
    starts <- seq_len(floor(lengths[[chrom]] / bin_size)) * bin_size - bin_size + 1
    GenomicRanges::GRanges(chrom, IRanges::IRanges(starts, width = bin_size))
  }))
  # Discard any bin overlapping a blacklist interval, not just endpoints.
  bins <- bins[!IRanges::overlapsAny(bins, blacklist, ignore.strand = TRUE)]
  if (!length(bins)) fail("No eligible nuclear bins after blacklisting")
  keys <- paste(as.character(GenomicRanges::seqnames(bins)),
                (IRanges::start(bins) - 1L) %/% bin_size, sep = ":")
  scores <- numeric(length(bins))
  connection <- gzfile(fragments, "rt"); on.exit(close(connection), add = TRUE)
  repeat {
    # Stream input to avoid loading a large fragments file into RAM.
    lines <- readLines(connection, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines)) next
    tab <- read.table(text = lines, sep = "\t", header = FALSE, quote = "",
      comment.char = "", colClasses = "character")
    if (ncol(tab) < 5L) fail("Fragments file requires at least five columns")
    tab <- tab[tab[[4]] %in% cells & tab[[1]] %in% chroms, , drop = FALSE]
    if (!nrow(tab)) next
    left <- suppressWarnings(as.numeric(tab[[2]])); right <- suppressWarnings(as.numeric(tab[[3]]))
    if (anyNA(c(left, right)) || any(left < 0 | right <= left) ||
        any(left != floor(left) | right != floor(right)) || any(right > lengths[tab[[1]]]))
      fail("Invalid fragment coordinates")
    # Each unique fragment contributes two endpoints; readSupport is deliberately
    # ignored because it contains PCR duplicates. BED end is exclusive: use end-1.
    hits <- match(c(paste(tab[[1]], left %/% bin_size, sep = ":"),
                    paste(tab[[1]], (right - 1) %/% bin_size, sep = ":")), keys)
    scores <- scores + tabulate(hits[!is.na(hits)], nbins = length(bins))
  }
  cutoff <- stats::quantile(scores, quantile_cutoff, names = FALSE)
  # Strict '>' follows the notebook. Ties can mean fewer than 1.5% are selected.
  selected <- bins[scores > cutoff]
  if (!length(selected)) fail("No bins above accessibility quantile; adjust --quantile")
  GenomicRanges::reduce(selected, ignore.strand = TRUE)
}

# Count sequenced A/C/G/T bases, NOT inferred insert lengths. CIGAR alignment
# blocks exclude deletions, introns, insertions and soft clips. Each base must
# pass the configured quality threshold and overlap the exact reference BED.
# Overlapping mates are counted separately; choose compatible mgatk settings.
count_read_bases <- function(cigar, position, sequence, quality, intervals, min_baseq) {
  # Newer Bioconductor versions moved CIGAR helpers to cigarillo. Retain
  # compatibility with the older GenomicAlignments in the exported environment.
  if (requireNamespace("cigarillo", quietly = TRUE)) {
    ref <- cigarillo::cigars_as_ranges_along_ref(cigar, lmmpos = position,
                                               ops = c("M", "=", "X"))[[1]]
    qry <- cigarillo::cigars_as_ranges_along_query(cigar, ops = c("M", "=", "X"))[[1]]
  } else {
    ref <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar, pos = position,
                                                           ops = c("M", "=", "X"))[[1]]
    qry <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar, ops = c("M", "=", "X"))[[1]]
  }
  if (!length(ref)) return(0)
  if (quality == "*" || !nzchar(quality)) return(0)
  bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
  qualities <- utf8ToInt(quality) - 33L
  if (length(qualities) != length(bases)) fail("Inconsistent BAM sequence/qualities")
  total <- 0
  for (block in seq_along(ref)) {
    rp <- seq.int(IRanges::start(ref)[block], IRanges::end(ref)[block])
    qp <- seq.int(IRanges::start(qry)[block], IRanges::end(qry)[block])
    # Reference intervals are reduced/sorted. Binary search finds the preceding
    # interval start, then its end determines whether each position is inside.
    ix <- findInterval(rp, IRanges::start(intervals))
    inside <- ix > 0L
    inside[inside] <- rp[inside] <= IRanges::end(intervals)[ix[inside]]
    total <- total + sum(inside & qualities[qp] >= min_baseq & bases[qp] %in% c("A", "C", "G", "T"))
  }
  total
}

nuclear_coverage <- function(path, reference, cells, barcode_lookup, mapq, baseq, keep_duplicates) {
  totals <- setNames(numeric(length(cells)), cells)
  reads <- totals
  ranges <- split(IRanges::ranges(reference), as.character(GenomicRanges::seqnames(reference)))
  bam <- Rsamtools::BamFile(path, yieldSize = 100000L)
  open(bam); on.exit(close(bam), add = TRUE)
  # Scan once in bounded chunks. Fetching each selected interval separately can
  # count one read twice when it overlaps multiple reference intervals.
  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
      isDuplicate = if (keep_duplicates) NA else FALSE),
    mapqFilter = mapq, tag = "CB", what = c("rname", "pos", "cigar", "seq", "qual"))
  repeat {
    chunk <- Rsamtools::scanBam(bam, param = param)[[1]]
    if (!length(chunk$pos)) break
    cb <- as.character(chunk$tag$CB)
    canonical <- unname(barcode_lookup[cb])
    chrom <- as.character(chunk$rname)
    valid <- which(!is.na(canonical) & canonical %in% cells & chrom %in% names(ranges))
    sequences <- as.character(chunk$seq); qualities <- as.character(chunk$qual)
    for (i in valid) {
      n <- count_read_bases(chunk$cigar[i], chunk$pos[i], sequences[i], qualities[i],
                            ranges[[chrom[i]]], baseq)
      if (n > 0) {
        totals[[canonical[i]]] <- totals[[canonical[i]]] + n
        reads[[canonical[i]]] <- reads[[canonical[i]]] + 1
      }
    }
    message("Scanned BAM chunk; accumulated nuclear bases: ", sum(totals))
  }
  list(bases = totals, reads = reads)
}

main <- function(args) {
  help <- paste(
    "Usage: Rscript scripts/11_estimate_mtdna_copy_number.R --sample ID",
    "  --bam ARC/outs/atac_possorted_bam.bam --barcodes retained_barcodes.txt",
    "  --mgatk-dir MGATK/final --out-prefix results/ID",
    "  (--reference-bed shared_reference.bed OR",
    "   --fragments ARC/outs/atac_fragments.tsv.gz --blacklist mm10-blacklist.bed)",
    "Optional: --barcode-map ARC/outs/per_barcode_metrics.csv",
    "  --bin-size 500 --quantile 0.985 --mapq 30 --baseq 20",
    "  --ploidy 2 --min-nuclear-bases 1000 --min-mt-depth 10",
    "  --keep-duplicates true|false (default false)",
    "Barcode file: one original ARC/Seurat barcode per line, no header.",
    "Match quality/duplicate conventions to the mgatk run; this script cannot infer them.", sep = "\n")
  if (!length(args) || "--help" %in% args) { cat(help, "\n"); return(invisible(NULL)) }
  if (length(args) %% 2L || any(!startsWith(args[seq(1L, length(args), 2L)], "--"))) fail(help)
  keys <- sub("^--", "", args[seq(1L, length(args), 2L)])
  if (anyDuplicated(keys)) fail("Repeated argument")
  opts <- as.list(args[seq(2L, length(args), 2L)]); names(opts) <- keys
  defaults <- list("bin-size" = "500", quantile = "0.985", mapq = "30", baseq = "20",
    ploidy = "2", "min-nuclear-bases" = "1000", "min-mt-depth" = "10", "keep-duplicates" = "false")
  allowed <- c(names(defaults), "sample", "bam", "barcodes", "mgatk-dir", "out-prefix",
               "reference-bed", "fragments", "blacklist", "barcode-map")
  if (any(!keys %in% allowed)) fail("Unknown argument: ", paste(setdiff(keys, allowed), collapse = ", "))
  for (key in names(defaults)) if (is.null(opts[[key]])) opts[[key]] <- defaults[[key]]
  for (key in c("sample", "bam", "barcodes", "mgatk-dir", "out-prefix"))
    if (is.null(opts[[key]]) || !nzchar(opts[[key]])) fail("Required: --", key)
  for (key in setdiff(names(defaults), "keep-duplicates")) {
    opts[[key]] <- suppressWarnings(as.numeric(opts[[key]]))
    if (length(opts[[key]]) != 1L || !is.finite(opts[[key]]) || opts[[key]] < 0) fail("Invalid --", key)
  }
  for (key in c("bin-size", "mapq", "baseq", "min-nuclear-bases"))
    if (opts[[key]] != floor(opts[[key]])) fail("--", key, " must be an integer")
  if (opts$quantile >= 1 || opts[["bin-size"]] < 1 || opts$ploidy <= 0) fail("Invalid quantile/bin size/ploidy")
  if (!opts[["keep-duplicates"]] %in% c("true", "false")) fail("Invalid --keep-duplicates")
  for (pkg in c("Rsamtools", "GenomicRanges", "GenomicAlignments"))
    if (!requireNamespace(pkg, quietly = TRUE)) fail("Missing R package: ", pkg)
  cells <- readLines(opts$barcodes, warn = FALSE)
  if (!length(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) fail("Empty/duplicate barcode list")
  lengths <- Rsamtools::scanBamHeader(opts$bam)[[1]]$targets
  reference <- if (!is.null(opts[["reference-bed"]])) {
    read_bed(opts[["reference-bed"]], lengths)
  } else {
    if (is.null(opts$fragments) || is.null(opts$blacklist)) fail("Provide reference BED or fragments AND blacklist")
    select_reference(opts$fragments, cells, lengths, read_bed(opts$blacklist, lengths),
                     opts[["bin-size"]], opts$quantile)
  }
  if (any(!grepl("^chr([1-9]|1[0-9])$", as.character(GenomicRanges::seqnames(reference)))))
    fail("Reference BED must contain mouse autosomes chr1-chr19 only")
  reference <- GenomicRanges::sort(GenomicRanges::reduce(reference, ignore.strand = TRUE))
  reference_bp <- sum(IRanges::width(reference))
  # Direct CB matching is preferred. If CB contains ATAC-library barcodes, the
  # ARC metrics table supplies the explicit ATAC->canonical barcode translation.
  lookup <- setNames(cells, cells)
  if (!is.null(opts[["barcode-map"]])) {
    map <- read.csv(opts[["barcode-map"]], colClasses = "character", check.names = FALSE)
    if (!all(c("barcode", "atac_barcode") %in% names(map))) fail("Barcode map needs barcode and atac_barcode")
    map <- map[map$barcode %in% cells, , drop = FALSE]
    if (anyNA(map[, c("barcode", "atac_barcode")]) || anyDuplicated(map$atac_barcode)) fail("Invalid barcode map")
    for (i in seq_len(nrow(map))) {
      old <- lookup[map$atac_barcode[i]]
      if (!is.na(old) && old != map$barcode[i]) fail("Conflicting barcode translation")
      lookup[map$atac_barcode[i]] <- map$barcode[i]
    }
  }
  # Use the genome-wide depth table, not coverage at m.5019. Missing mgatk cells
  # stay NA: lack of a measurement must never be interpreted as zero copies.
  depth_files <- list.files(opts[["mgatk-dir"]], pattern = "[.]depthTable[.]txt$", full.names = TRUE)
  if (length(depth_files) != 1L) fail("Expected exactly one mgatk depthTable.txt")
  depth <- read.table(depth_files, header = FALSE, colClasses = c("character", "numeric"))
  if (ncol(depth) != 2L || anyDuplicated(depth[[1]]) || anyNA(depth) ||
      any(!is.finite(depth[[2]]) | depth[[2]] < 0)) fail("Invalid mgatk depth table")
  # Apply the same explicit translation if mgatk was run with ATAC barcodes.
  canonical_mt <- unname(lookup[depth[[1]]])
  matched <- !is.na(canonical_mt) & canonical_mt %in% cells
  if (anyDuplicated(canonical_mt[matched])) fail("Multiple mgatk rows map to the same retained cell")
  mt_depth <- depth[[2]][matched][match(cells, canonical_mt[matched])]
  if (all(is.na(mt_depth))) fail("No retained barcodes match mgatk depth table")
  measured <- nuclear_coverage(opts$bam, reference, cells, lookup, opts$mapq, opts$baseq,
                               opts[["keep-duplicates"]] == "true")
  if (!sum(measured$bases)) fail("No nuclear bases matched: check CB tags, mapping and reference BED")
  nuclear_depth <- unname(measured$bases) / reference_bp
  estimate <- ifelse(nuclear_depth > 0, opts$ploidy * mt_depth / nuclear_depth, NA_real_)
  status <- ifelse(is.na(mt_depth), "missing_mgatk",
    ifelse(measured$bases == 0, "zero_nuclear_coverage",
      ifelse(measured$bases < opts[["min-nuclear-bases"]], "low_nuclear_coverage",
        ifelse(mt_depth < opts[["min-mt-depth"]], "low_mitochondrial_coverage", "pass"))))
  out <- data.frame(sample_id = opts$sample, barcode = cells,
    nuclear_reference_bp = reference_bp, nuclear_aligned_bases = unname(measured$bases),
    nuclear_reads_with_counted_bases = unname(measured$reads), nuclear_mean_depth = nuclear_depth,
    mitochondrial_mean_depth = mt_depth, assumed_nuclear_ploidy = opts$ploidy,
    mtDNA_copy_estimate_raw = estimate,
    mtDNA_copy_estimate = ifelse(status == "pass", estimate, NA_real_), qc_status = status)
  prefix <- opts[["out-prefix"]]
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  write.table(out, paste0(prefix, ".mtDNA_copy_number.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  bed <- data.frame(chrom = as.character(GenomicRanges::seqnames(reference)),
                    start = IRanges::start(reference) - 1L, end = IRanges::end(reference))
  write.table(bed, paste0(prefix, ".nuclear_reference.bed"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  # Save all thresholds and versions alongside results for reproducibility.
  capture.output(list(arguments = opts, nuclear_reference_bp = reference_bp,
    mgatk_depth_file = depth_files, overlap_policy = "count each mate separately",
    session = sessionInfo()), file = paste0(prefix, ".settings.txt"))
  message("Wrote estimates: ", prefix, ".mtDNA_copy_number.tsv")
  invisible(out)
}

# Sourcing exposes helpers for testing without executing the command-line entry.
if (sys.nframe() == 0L) {
  tryCatch(main(commandArgs(trailingOnly = TRUE)), error = function(e) {
    message("ERROR: ", conditionMessage(e)); quit(status = 1L)
  })
}
