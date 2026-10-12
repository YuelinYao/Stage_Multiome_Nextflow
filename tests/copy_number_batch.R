# Compare the batched counter against the original scalar implementation.
# Run: Rscript tests/copy_number_batch.R [--benchmark]
source('scripts/10_estimate_mtdna_copy_number.R')
set.seed(82)
cigars <- c('40M', '5S20M3I10M2S', '10M5D15M', '6M100N6M',
            '5H10M5H', '10=5X15M', '5S', '20M')
query_length <- function(cigar) {
  tokens <- regmatches(cigar, gregexpr('[0-9]+[MIDNSHP=X]', cigar))[[1]]
  sum(as.numeric(sub('[A-Z=]$', '', tokens[grepl('[MIS=X]$', tokens)])))
}
reference <- GenomicRanges::GRanges(c('chr1', 'chr1', 'chr1', 'chr2'),
  IRanges::IRanges(c(1, 20, 120, 60), c(8, 35, 170, 110)))
ranges <- split(IRanges::ranges(reference), as.character(GenomicRanges::seqnames(reference)))
make_reads <- function(n) {
  cigar <- sample(cigars, n, replace = TRUE)
  lengths <- vapply(cigar, query_length, numeric(1))
  sequence <- vapply(lengths, function(n) paste(sample(c('A', 'C', 'G', 'T', 'N'),
    n, replace = TRUE), collapse = ''), character(1))
  quality <- vapply(lengths, function(n) paste(sample(c('I', '5', '!', '?', ' '),
    n, replace = TRUE, prob = c(.85, .05, .04, .05, .01)), collapse = ''), character(1))
  # Missing BAM quality values must retain zero coverage.
  quality[seq_len(min(5, n))] <- '*'
  list(cigar = cigar, position = sample(1:190, n, replace = TRUE),
       chromosome = sample(c('chr1', 'chr2'), n, replace = TRUE),
       sequence = sequence, quality = quality)
}
scalar <- function(x, baseq) vapply(seq_along(x$cigar), function(i)
  count_read_bases(x$cigar[i], x$position[i], x$sequence[i], x$quality[i],
                   ranges[[x$chromosome[i]]], baseq), numeric(1))
batch <- function(x, baseq) count_alignment_bases(x$cigar, x$position, x$chromosome,
  Biostrings::DNAStringSet(x$sequence), Biostrings::BStringSet(x$quality), reference, baseq)
x <- make_reads(1000)
for (baseq in c(0, 9, 20, 40, 94, 100))
  stopifnot(identical(scalar(x, baseq), batch(x, baseq)))
# CIGAR deletion and skipped-reference blocks contribute no bases, even when
# they cross selected BED intervals. Empty batches remain valid.
stopifnot(identical(batch(make_reads(0), 20), numeric()))
cat('Batched and scalar counts agree for mixed CIGARs, qualities and reference intervals.\n')

# Cross a real 100,000-record BAM chunk boundary. An initial chunk with no
# retained cells must not prevent later chunks from contributing coverage.
check_bam_chunks <- function() {
  root <- tempfile('copy-batch-bam-'); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  sam <- file.path(root, 'chunks.sam')
  writeLines(c('@HD\tVN:1.6\tSO:coordinate', '@SQ\tSN:chr1\tLN:200',
    sprintf('other%d\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tCB:Z:OTHER',
            seq_len(100000)),
    'one\t0\tchr1\t50\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tCB:Z:CELL1',
    'missing\t0\tchr1\t50\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\tCB:Z:CELL1',
    'two\t0\tchr1\t50\t60\t5M2I5M\t*\t0\t0\tAAAAAAAAAAAA\tIIIIIIIIIIII\tCB:Z:CELL2'), sam)
  bam <- Rsamtools::asBam(sam, destination = file.path(root, 'chunks'), overwrite = TRUE)
  ref <- GenomicRanges::GRanges('chr1', IRanges::IRanges(c(50, 57), c(52, 59)))
  cells <- c('CELL1', 'CELL2')
  for (baseq in c(0, 20)) {
    measured <- nuclear_coverage(bam, ref, cells, setNames(cells, cells), 30, baseq, FALSE)
    stopifnot(identical(unname(measured$bases), c(6, 6)),
              identical(unname(measured$reads), c(1, 1)))
  }
}
check_bam_chunks()
cat('BAM chunk boundaries, disjoint intervals and missing qualities passed.\n')

if ('--benchmark' %in% commandArgs(trailingOnly = TRUE)) {
  x <- make_reads(20000)
  old <- system.time(expected <- scalar(x, 20))['elapsed']
  new <- system.time(actual <- batch(x, 20))['elapsed']
  stopifnot(identical(expected, actual))
  cat(sprintf('20,000 synthetic alignments: scalar %.3f s; batch %.3f s; %.1fx faster.\n',
              old, new, old / new))
  cat('This benchmarks counting only, excluding BAM/fragments I/O and reference selection.\n')
}
