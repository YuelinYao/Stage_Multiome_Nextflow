# Run from the repository root; no Signac or real sequencing inputs required.
for (f in list.files("scripts", pattern = "[.]R$", full.names = TRUE)) parse(f)
env <- new.env(parent = globalenv())
# Load pure helpers without loading the unavailable scientific dependencies.
for (expr in parse("scripts/00_common.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-"))) eval(expr, env)
}
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
cfg <- jsonlite::fromJSON("examples/qc_thresholds.example.json", simplifyVector = FALSE)
config_file <- tempfile(fileext = ".json")
write_config <- function(x) {
  jsonlite::write_json(x, config_file, auto_unbox = TRUE, null = "null")
  config_file
}
valid <- env$read_analysis_config(write_config(cfg))
stopifnot(identical(valid$clustering$RNA_dims, 1:20))
bad <- cfg; bad$sample_overrides <- list(list(nCount_RNA = list(min = 1, max = 2)))
expect_error(env$read_analysis_config(write_config(bad)))
bad <- cfg; bad$thresholds$nCount_RNA <- list(min = 2, max = 1)
expect_error(env$read_analysis_config(write_config(bad)))
bad <- cfg; bad$clustering$RNA_dims <- "20:1"
expect_error(env$read_analysis_config(write_config(bad)))
tab <- data.frame(barcode = c("lower", "upper", "below", "missing", "infinite"))
limits <- setNames(lapply(env$qc_metric_names, function(x) list(min = 1, max = 2)),
                   env$qc_metric_names)
for (metric in env$qc_metric_names) tab[[metric]] <- c(1, 2, 0, NA, Inf)
stopifnot(identical(env$apply_qc(tab, limits), c("lower", "upper")))

# Execute the mgatk extractor with a mock ReadMGATK and a plain input object.
obj_file <- tempfile(); out_file <- tempfile()
obj <- matrix(0, 1, 2, dimnames = list("gene", c("sample_present", "sample_missing")))
class(obj) <- c("mock_cells", class(obj))
`$.mock_cells` <- function(x, name) c("present", "missing")
saveRDS(obj, obj_file)
env$ReadMGATK <- function(dir) list(
  refallele = data.frame(pos = 5019, ref = "A"),
  counts = matrix(c(3, 3, 2, 2, 0, 0, 0, 0), ncol = 1,
    dimnames = list(paste(rep(c("A", "G", "C", "T"), each = 2), 5019,
                          rep(c("fwd", "rev"), 4), sep = "-"), "present")),
  depth = data.frame(depth = 100, row.names = "present"))
env$args <- c("sample", obj_file, "mock", out_file)
exprs <- parse("scripts/06c_mgatk_metrics.R")
for (expr in exprs[-(1:2)]) eval(expr, env)
result <- read.delim(out_file)
stopifnot(result$hf_5019[1] == 0.4, result$depth_5019[1] == 10,
          result$mtDNA_depth[1] == 100, is.na(result$hf_5019[2]),
          is.na(result$mtDNA_depth[2]))
unlink(c(config_file, obj_file, out_file))
cat("All regression checks passed.\n")
