#!/usr/bin/env Rscript
# Estimate mtDNA abundance relative to accessible, diploid nuclear DNA.
#
# Formula: estimated copies = nuclear ploidy * mean mtDNA depth / mean nuclear depth.
# This is an assay-dependent estimate, not an independently calibrated molecule count.
# main.nf runs this per sample from step 06 singlets, alongside step 09.
# It can also run directly.
# Run with --help for arguments. All BED files use 0-based, end-exclusive coordinates.
# Required R packages are listed in environment.yml. Signac is needed only
# when selecting bins with its default blacklist; a custom BED avoids that dependency.

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

# Load the bundled GRanges directly: no download or intermediate BED is needed.
# Keep a custom --blacklist option because newer Signac releases may remove
# this deprecated dataset. A supplied reference BED bypasses blacklist loading.
load_signac_mm10_blacklist <- function() {
  if (!requireNamespace("Signac", quietly = TRUE))
    fail("Install Signac for its mm10 blacklist, or supply --blacklist /path/to/mm10-blacklist.bed")
  blacklist <- tryCatch(getExportedValue("Signac", "blacklist_mm10"),
    error = function(e) fail("Signac::blacklist_mm10 is unavailable; supply --blacklist /path/to/mm10-blacklist.bed"))
  if (!inherits(blacklist, "GRanges") || !length(blacklist))
    fail("Signac::blacklist_mm10 must contain genomic ranges")
  blacklist
}

# Score only full-width autosomal bins. Sex chromosomes are excluded so a
# diploid normalization does not depend on inferred sex. The chromosome lengths
# come from YOUR BAM header, rather than a hard-coded human or mouse genome.
select_reference <- function(fragments, cells, lengths, blacklist, bin_size, quantile_cutoff) {
  chroms <- grep("^chr([1-9]|1[0-9])$", names(lengths), value = TRUE)
  if (!length(chroms)) fail("Reference selection requires mouse chr1-chr19 contigs in the BAM header")
  bin_counts <- floor(lengths[chroms] / bin_size)
  chroms <- chroms[bin_counts > 0]
  bin_counts <- bin_counts[chroms]
  if (!length(chroms)) fail("No full-width autosomal bins; reduce --bin-size")
  offsets <- setNames(c(0, head(cumsum(bin_counts), -1L)), chroms)
  bins <- do.call(c, lapply(chroms, function(chrom) {
    starts <- seq_len(floor(lengths[[chrom]] / bin_size)) * bin_size - bin_size + 1
    GenomicRanges::GRanges(chrom, IRanges::IRanges(starts, width = bin_size))
  }))
  # Discard any bin overlapping a blacklist interval, not just endpoints.
  if (length(blacklist))
    bins <- bins[!IRanges::overlapsAny(bins, blacklist, ignore.strand = TRUE)]
  if (!length(bins)) fail("No eligible nuclear bins after blacklisting")
  # Integer lookup avoids building and repeatedly hashing millions of strings
  # such as 'chr1:1234'. Zero marks a blacklisted (ineligible) bin.
  bin_lookup <- integer(sum(bin_counts))
  index <- unname(offsets[as.character(GenomicRanges::seqnames(bins))]) +
    (IRanges::start(bins) - 1L) %/% bin_size + 1L
  bin_lookup[index] <- seq_along(bins)
  scores <- numeric(length(bins))
  message("Selecting nuclear reference from ", length(bins), " eligible bins")
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
    local <- c(left %/% bin_size, (right - 1) %/% bin_size)
    # A chromosome's trailing partial bin is intentionally excluded. Check its
    # local index before applying offsets, so it cannot alias the next chromosome.
    full_width <- local < rep(unname(bin_counts[tab[[1]]]), 2L)
    global <- local[full_width] + rep(unname(offsets[tab[[1]]]), 2L)[full_width] + 1L
    scores <- scores + tabulate(bin_lookup[global], nbins = length(bins))
  }
  cutoff <- stats::quantile(scores, quantile_cutoff, names = FALSE)
  # Strict '>' follows the notebook. Ties can mean fewer than 1.5% are selected.
  selected <- bins[scores > cutoff]
  if (!length(selected)) fail("No bins above accessibility quantile; adjust --quantile")
  message("Selected ", length(selected), " bins above score cutoff ", cutoff)
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
  if (is.na(quality) || quality == "*" || !nzchar(quality) ||
      is.na(sequence) || sequence == "*" || !nzchar(sequence)) return(0)
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

# Return one count per alignment. Parse CIGARs and intersect reference blocks
# in batches, rather than creating a vector of every aligned base for every read.
# Only M/= /X blocks contribute; query offsets preserve insertions and clips.
count_alignment_bases <- function(cigar, position, chromosome, sequence, quality,
                                  reference, min_baseq,
                                  reference_index = GenomicRanges::GNCList(reference)) {
  counts <- numeric(length(cigar))
  if (!length(cigar)) return(counts)
  if (requireNamespace("cigarillo", quietly = TRUE)) {
    ref <- cigarillo::cigars_as_ranges_along_ref(cigar, lmmpos = position,
                                               ops = c("M", "=", "X"))
    qry <- cigarillo::cigars_as_ranges_along_query(cigar, ops = c("M", "=", "X"))
  } else {
    ref <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar, pos = position,
                                                           ops = c("M", "=", "X"))
    qry <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar, ops = c("M", "=", "X"))
  }
  block_read <- rep.int(seq_along(cigar), lengths(ref))
  if (!length(block_read)) return(counts)
  ref <- unlist(ref, use.names = FALSE)
  qry <- unlist(qry, use.names = FALSE)
  blocks <- GenomicRanges::GRanges(chromosome[block_read], ref)
  hits <- GenomicRanges::findOverlaps(blocks, reference_index, ignore.strand = TRUE)
  if (!length(hits)) return(counts)
  block <- S4Vectors::queryHits(hits)
  region <- S4Vectors::subjectHits(hits)
  read <- block_read[block]

  # Intersections give the exact query subsequences inside the BED. One read
  # can contribute several pieces, but will still be counted as one read below.
  left <- pmax(IRanges::start(ref)[block], IRanges::start(reference)[region])
  right <- pmin(IRanges::end(ref)[block], IRanges::end(reference)[region])
  query_start <- IRanges::start(qry)[block] + left - IRanges::start(ref)[block]
  query_end <- query_start + right - left
  sequence <- sequence[read]
  quality <- methods::as(quality[read], "BStringSet")
  usable <- Biostrings::width(sequence) > 0L & Biostrings::width(quality) > 0L
  # Preserve the original missing-quality rule, including the SAM '*' sentinel.
  single <- which(Biostrings::width(quality) == 1L)
  usable[single] <- usable[single] & as.character(quality[single]) != "*"
  if (!any(usable)) return(counts)
  if (any(Biostrings::width(sequence)[usable] != Biostrings::width(quality)[usable]))
    fail("Inconsistent BAM sequence/qualities")
  read <- read[usable]
  sequence <- Biostrings::subseq(sequence[usable], start = query_start[usable],
                                end = query_end[usable])
  quality <- Biostrings::subseq(quality[usable], start = query_start[usable],
                               end = query_end[usable])
  piece_counts <- as.numeric(Biostrings::letterFrequency(sequence, letters = "ACGT"))
  if (min_baseq >= 0) {
    # PHRED+33: include the negative-quality encoding too. Rsamtools converts
    # absent BAM qualities to spaces (score -1), which must fail even baseq=0.
    low_letters <- paste(intToUtf8(seq.int(1L, min(126L, 32L + min_baseq)),
                                  multiple = TRUE), collapse = "")
    low_quality <- which(Biostrings::letterFrequency(quality, letters = low_letters)[, 1] > 0)
    # Most pieces can use compiled letterFrequency directly. For pieces with
    # low-quality bases, intersect the base and quality masks at each position.
    sequences <- as.character(sequence[low_quality])
    qualities <- as.character(quality[low_quality])
    for (i in seq_along(low_quality)) {
      bases <- strsplit(sequences[i], "", fixed = TRUE)[[1]]
      piece_counts[low_quality[i]] <- sum(bases %in% c("A", "C", "G", "T") &
                                          utf8ToInt(qualities[i]) - 33L >= min_baseq)
    }
  }
  summed <- rowsum(piece_counts, group = read, reorder = FALSE)
  counts[as.integer(rownames(summed))] <- summed[, 1]
  counts
}

nuclear_coverage <- function(path, reference, cells, barcode_lookup, mapq, baseq, keep_duplicates) {
  totals <- setNames(numeric(length(cells)), cells)
  reads <- totals
  reference <- GenomicRanges::sort(GenomicRanges::reduce(reference, ignore.strand = TRUE))
  reference_index <- GenomicRanges::GNCList(reference)
  chromosomes <- unique(as.character(GenomicRanges::seqnames(reference)))
  bam <- Rsamtools::BamFile(path, yieldSize = 100000L)
  open(bam); on.exit(close(bam), add = TRUE)
  # Scan once in bounded chunks. Fetching each selected interval separately can
  # count one read twice when it overlaps multiple reference intervals.
  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
      isDuplicate = if (keep_duplicates) NA else FALSE),
    mapqFilter = mapq, tag = "CB", what = c("rname", "pos", "cigar", "seq", "qual"))
  chunk_number <- 0L
  repeat {
    chunk <- Rsamtools::scanBam(bam, param = param)[[1]]
    if (!length(chunk$pos)) break
    chunk_number <- chunk_number + 1L
    cb <- as.character(chunk$tag$CB)
    canonical <- unname(barcode_lookup[cb])
    chrom <- as.character(chunk$rname)
    cell_index <- match(canonical, cells)
    valid <- which(!is.na(cell_index) & chrom %in% chromosomes)
    if (length(valid)) {
      n <- count_alignment_bases(chunk$cigar[valid], chunk$pos[valid], chrom[valid],
        chunk$seq[valid], chunk$qual[valid], reference, baseq, reference_index)
      grouped <- rowsum(n, group = cell_index[valid], reorder = FALSE)
      index <- as.integer(rownames(grouped))
      totals[index] <- totals[index] + grouped[, 1]
      # A read spanning several selected regions contributes one read count.
      reads <- reads + tabulate(cell_index[valid][n > 0], nbins = length(cells))
    }
    message("Scanned BAM chunk ", chunk_number, ": ", length(chunk$pos),
      " filtered alignments; ", length(valid),
      " retained-cell nuclear alignments; accumulated nuclear bases: ", sum(totals))
  }
  list(bases = totals, reads = reads)
}

# Each source is aligned independently to the retained canonical barcodes.
# Depth is already a genome-wide mean: do not divide by mtDNA length again.
read_mgatk_depth <- function(directory, cells, lookup) {
  if (!dir.exists(directory)) fail("Missing mgatk directory: ", directory)
  files <- list.files(directory, pattern = "[.]depthTable[.]txt$", full.names = TRUE)
  if (length(files) != 1L) fail("Expected exactly one mgatk depthTable.txt in: ", directory)
  depth <- read.table(files, header = FALSE, colClasses = c("character", "numeric"))
  if (ncol(depth) != 2L || anyDuplicated(depth[[1]]) || anyNA(depth) ||
      any(!nzchar(depth[[1]])) || any(!is.finite(depth[[2]]) | depth[[2]] < 0))
    fail("Invalid mgatk depth table: ", files)
  canonical <- unname(lookup[depth[[1]]])
  matched <- !is.na(canonical) & canonical %in% cells
  if (anyDuplicated(canonical[matched])) fail("Multiple mgatk rows map to the same retained cell: ", files)
  # Unmatched cells remain NA, including when this entire source lacks matches.
  # This must not suppress a valid estimate from the other source.
  list(depth = depth[[2]][matched][match(cells, canonical[matched])], file = files)
}

main <- function(args) {
  help <- paste(
    "Usage: Rscript scripts/10_estimate_mtdna_copy_number.R --sample ID",
    "  --bam ARC/outs/atac_possorted_bam.bam",
    "  (--barcodes retained_barcodes.txt OR --metrics-tsv sample.mtDNA.tsv)",
    "  --mgatk-dir MGATK/final --mgatk-enriched-dir MGATK_ENRICHED/final",
    "  --out-prefix results/ID (provide either or both mgatk directories)",
    "  (--reference-bed shared_reference.bed OR",
    "   --fragments ARC/outs/atac_fragments.tsv.gz)",
    "Optional: --json pipeline_settings.json (copy_number settings)",
    "  --barcode-map ARC/outs/per_barcode_metrics.csv",
    "  --blacklist mm10-blacklist.bed (default: Signac::blacklist_mm10)",
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
  allowed <- c(names(defaults), "sample", "bam", "barcodes", "mgatk-dir", "mgatk-enriched-dir", "out-prefix",
               "reference-bed", "fragments", "blacklist", "barcode-map", "metrics-tsv", "json")
  if (any(!keys %in% allowed)) fail("Unknown argument: ", paste(setdiff(keys, allowed), collapse = ", "))
  # JSON supplies analysis settings; explicit script arguments take precedence.
  if (!is.null(opts$json)) {
    settings <- jsonlite::fromJSON(opts$json, simplifyVector = FALSE)
    supplied <- settings$copy_number
    if (!is.null(supplied)) {
      fields <- c("bin_size", "quantile", "mapq", "baseq", "ploidy",
        "min_nuclear_bases", "min_mt_depth", "keep_duplicates", "reference_bed", "blacklist")
      if (!is.list(supplied) || (length(supplied) &&
          (is.null(names(supplied)) || anyDuplicated(names(supplied)))) ||
          any(!names(supplied) %in% fields)) fail("Invalid copy_number JSON settings")
      for (field in names(supplied)) {
        value <- supplied[[field]]
        if (is.null(value)) next
        key <- gsub("_", "-", field, fixed = TRUE)
        if (field == "keep_duplicates") {
          if (length(value) != 1L || !is.logical(value) || is.na(value)) fail("Invalid copy_number.", field)
          value <- if (value) "true" else "false"
        } else if (field %in% c("reference_bed", "blacklist")) {
          if (length(value) != 1L || !is.character(value) || is.na(value) || !startsWith(value, "/"))
            fail("Invalid copy_number.", field)
        } else if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
          fail("Invalid copy_number.", field)
        }
        if (is.null(opts[[key]])) opts[[key]] <- as.character(value)
      }
    }
  }
  for (key in names(defaults)) if (is.null(opts[[key]])) opts[[key]] <- defaults[[key]]
  for (key in c("sample", "bam", "out-prefix"))
    if (is.null(opts[[key]]) || !nzchar(opts[[key]])) fail("Required: --", key)
  source_keys <- c(mgatk = "mgatk-dir", enriched = "mgatk-enriched-dir")
  source_keys <- source_keys[vapply(source_keys, function(key) !is.null(opts[[key]]), logical(1))]
  if (!length(source_keys)) fail("Provide --mgatk-dir and/or --mgatk-enriched-dir")
  for (key in source_keys) if (!nzchar(opts[[key]])) fail("Empty --", key)
  for (key in setdiff(names(defaults), "keep-duplicates")) {
    opts[[key]] <- suppressWarnings(as.numeric(opts[[key]]))
    if (length(opts[[key]]) != 1L || !is.finite(opts[[key]]) || opts[[key]] < 0) fail("Invalid --", key)
  }
  for (key in c("bin-size", "mapq", "baseq", "min-nuclear-bases"))
    if (opts[[key]] != floor(opts[[key]])) fail("--", key, " must be an integer")
  if (opts$quantile >= 1 || opts[["bin-size"]] < 1 || opts$ploidy <= 0) fail("Invalid quantile/bin size/ploidy")
  if (!opts[["keep-duplicates"]] %in% c("true", "false")) fail("Invalid --keep-duplicates")
  for (pkg in c("Rsamtools", "GenomicRanges", "GenomicAlignments", "Biostrings"))
    if (!requireNamespace(pkg, quietly = TRUE)) fail("Missing R package: ", pkg)
  if (sum(!vapply(opts[c("barcodes", "metrics-tsv")], is.null, logical(1))) != 1L)
    fail("Provide exactly one of --barcodes and --metrics-tsv")
  cells <- if (!is.null(opts[["metrics-tsv"]])) {
    metrics <- read.delim(opts[["metrics-tsv"]], check.names = FALSE,
      colClasses = c(sample_id = "character", barcode = "character"))
    if (!all(c("sample_id", "barcode") %in% names(metrics)) || anyDuplicated(names(metrics)) ||
        !nrow(metrics) || anyNA(metrics$sample_id) || any(metrics$sample_id != opts$sample))
      fail("Invalid sample/columns in --metrics-tsv")
    metrics$barcode
  } else readLines(opts$barcodes, warn = FALSE)
  if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) fail("Empty/duplicate barcode list")
  lengths <- Rsamtools::scanBamHeader(opts$bam)[[1]]$targets
  # Direct CB matching is preferred. If CB contains ATAC-library barcodes, the
  # ARC metrics table supplies the explicit ATAC->canonical barcode translation.
  lookup <- setNames(cells, cells)
  if (!is.null(opts[["barcode-map"]])) {
    map <- read.csv(opts[["barcode-map"]], colClasses = "character", check.names = FALSE)
    if (!all(c("barcode", "atac_barcode") %in% names(map)) || anyDuplicated(names(map)))
      fail("Barcode map needs unique barcode and atac_barcode columns")
    map <- map[map$barcode %in% cells, , drop = FALSE]
    if (anyNA(map[, c("barcode", "atac_barcode")]) ||
        any(!nzchar(map$barcode)) || any(!nzchar(map$atac_barcode)) ||
        anyDuplicated(map$barcode) || anyDuplicated(map$atac_barcode))
      fail("Barcode map must have nonempty, one-to-one barcode mappings")
    for (i in seq_len(nrow(map))) {
      old <- lookup[map$atac_barcode[i]]
      if (!is.na(old) && old != map$barcode[i]) fail("Conflicting barcode translation")
      lookup[map$atac_barcode[i]] <- map$barcode[i]
    }
  }
  blacklist_source <- "not applied: supplied reference BED"
  reference <- if (!is.null(opts[["reference-bed"]])) {
    read_bed(opts[["reference-bed"]], lengths)
  } else {
    if (is.null(opts$fragments)) fail("Provide --reference-bed or --fragments")
    blacklist <- if (!is.null(opts$blacklist)) {
      blacklist_source <- normalizePath(opts$blacklist, mustWork = TRUE)
      read_bed(opts$blacklist, lengths)
    } else {
      blacklist_source <- "Signac::blacklist_mm10"
      bundled <- load_signac_mm10_blacklist()
      # The bundled dataset covers the complete genome, while the BAM may
      # contain a subset of its contigs. Unused contigs cannot overlap our bins.
      bundled[as.character(GenomicRanges::seqnames(bundled)) %in% names(lengths)]
    }
    # Fragment barcodes may be canonical ARC IDs or their ATAC aliases.
    # Apply the same accepted IDs used for BAM and mgatk depth matching.
    select_reference(opts$fragments, names(lookup), lengths, blacklist,
                     opts[["bin-size"]], opts$quantile)
  }
  if (any(!grepl("^chr([1-9]|1[0-9])$", as.character(GenomicRanges::seqnames(reference)))))
    fail("Reference BED must contain mouse autosomes chr1-chr19 only")
  reference <- GenomicRanges::sort(GenomicRanges::reduce(reference, ignore.strand = TRUE))
  reference_bp <- sum(IRanges::width(reference))
  # Read both mitochondrial sources before scanning the nuclear BAM once.
  # Both ratios use exactly the same retained cells and nuclear reference.
  depths <- lapply(source_keys, function(key) read_mgatk_depth(opts[[key]], cells, lookup))
  measured <- nuclear_coverage(opts$bam, reference, cells, lookup, opts$mapq, opts$baseq,
                               opts[["keep-duplicates"]] == "true")
  if (!sum(measured$bases)) fail("No nuclear bases matched: check CB tags, mapping and reference BED")
  nuclear_depth <- unname(measured$bases) / reference_bp
  out <- data.frame(sample_id = opts$sample, barcode = cells,
    nuclear_reference_bp = reference_bp, nuclear_aligned_bases = unname(measured$bases),
    nuclear_reads_with_counted_bases = unname(measured$reads), nuclear_mean_depth = nuclear_depth,
    assumed_nuclear_ploidy = opts$ploidy)
  # Always expose both requested metadata names, including when the ordinary
  # source has not been supplied yet. Missing sources must never become zero.
  for (suffix in c("", "_enriched")) {
    out[[paste0("mitochondrial_mean_depth", suffix)]] <- NA_real_
    out[[paste0("mtDNA_CN_normalised", suffix, "_raw")]] <- NA_real_
    out[[paste0("mtDNA_CN_normalised", suffix)]] <- NA_real_
    out[[paste0("qc_status", suffix)]] <- "missing_mgatk_source"
  }
  for (source in names(depths)) {
    mt_depth <- depths[[source]]$depth
    estimate <- ifelse(nuclear_depth > 0, opts$ploidy * mt_depth / nuclear_depth, NA_real_)
    status <- ifelse(is.na(mt_depth), "missing_mgatk",
      ifelse(measured$bases == 0, "zero_nuclear_coverage",
        ifelse(measured$bases < opts[["min-nuclear-bases"]], "low_nuclear_coverage",
          ifelse(mt_depth < opts[["min-mt-depth"]], "low_mitochondrial_coverage", "pass"))))
    # The unsuffixed estimate uses ordinary mgatk; _enriched uses enriched mgatk.
    suffix <- if (source == "enriched") "_enriched" else ""
    out[[paste0("mitochondrial_mean_depth", suffix)]] <- mt_depth
    out[[paste0("mtDNA_CN_normalised", suffix, "_raw")]] <- estimate
    out[[paste0("mtDNA_CN_normalised", suffix)]] <- ifelse(status == "pass", estimate, NA_real_)
    out[[paste0("qc_status", suffix)]] <- unname(status)
    # Preserve previous ordinary-mgatk column names for existing consumers.
    if (source == "mgatk") {
      out$mtDNA_copy_estimate_raw <- estimate
      out$mtDNA_copy_estimate <- out$mtDNA_CN_normalised
    }
  }
  prefix <- opts[["out-prefix"]]
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  write.table(out, paste0(prefix, ".mtDNA_copy_number.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  bed <- data.frame(chrom = as.character(GenomicRanges::seqnames(reference)),
                    start = IRanges::start(reference) - 1L, end = IRanges::end(reference))
  write.table(bed, paste0(prefix, ".nuclear_reference.bed"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  # Save all thresholds and versions alongside results for reproducibility.
  capture.output(list(arguments = opts, nuclear_reference_bp = reference_bp,
    mgatk_depth_files = lapply(depths, function(x) x$file), blacklist_source = blacklist_source,
    coverage_method = "batched CIGAR intersections and quality-filtered sequence counts",
    overlap_policy = "count each mate separately",
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
