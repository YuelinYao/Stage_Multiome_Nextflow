# Per-stage m.5019A>G Multiome workflow

This Nextflow workflow analyses one mouse developmental stage at a time (mm10). It calls peaks from each sample's ATAC fragments, creates one common peak set for the stage, reviews eight RNA/ATAC QC metrics, removes reviewed low-quality clusters and RNA doublets per sample, preprocesses and integrates RNA and ATAC, builds WNN clusters, attaches precomputed mgatk results, and identifies RNA markers for each WNN cluster.

## Requirements and inputs

Use Nextflow with Java 17 or later. Load Java before installing or running Nextflow:

```bash
# Example (Cambridge MBU cluster)
module add java/17
curl -s https://get.nextflow.io | bash
```

The installation command creates a `nextflow` executable in the current directory. Use `./nextflow` there, or place it in a directory on your `PATH` to use the `nextflow` commands below.

Use these existing environments when running this project:

| Command-line input | Path |
|---|---|
| `--conda_env` | `/suffolk/WorkGenomicsD/yy605/conda_env/stage_multiome` |
| `--macs3_path` | `/suffolk/WorkGenomicsD/yy605/conda_env/macs3/bin/macs3` |

Enable the Conda profile with `-profile conda` for local execution or `-profile slurm,conda` for Slurm. Nextflow uses the R environment above for each task and calls MACS3 through its explicit executable path. Both paths must be accessible on every execution node.

Package versions for the existing `stage_multiome` Conda environment are recorded in [environment.yml](environment.yml). This is an exported Linux environment, including its dependencies and installation prefix.

The sample CSV must contain only `stage,sample_id,arc_out,mgatk_dir`; see [examples/samples.csv](examples/samples.csv). For example:

```csv
stage,sample_id,arc_out,mgatk_dir
E8.5,5019_Hi_1,/path/to/arc/outs,/path/to/mgatk/final
```

Sample IDs must be unique within each stage, start with a letter or number, and contain only letters, numbers, `_`, `.` or `-`. Sample IDs are workflow labels and do not need to match upstream Cell Ranger or mgatk run names; the paths must identify the corresponding sample outputs, and sample-specific QC JSON keys must match these IDs. All execution nodes must be able to read the absolute input paths:

| Input | Required content |
|---|---|
| `arc_out` | Cell Ranger ARC `outs` directory containing `filtered_feature_bc_matrix.h5`, `atac_fragments.tsv.gz` and its `.tbi` index |
| `mgatk_dir` | Completed mgatk output readable by `Signac::ReadMGATK`, usually the `final` directory |

Gene symbols, the reference genome and barcodes are requirements for the existing data, not additional files or CSV columns:

| Data requirement | Where it comes from |
|---|---|
| Mouse gene symbols | Gene names stored in ARC's `filtered_feature_bc_matrix.h5`; cell-cycle scoring and sex inference require symbols rather than Ensembl IDs |
| mm10 reference | ARC and mgatk outputs must use a compatible mm10 reference; the workflow supplies mm10 gene annotation through `EnsDb.Mmusculus.v79` and checks that mgatk's reference allele at mitochondrial position 5019 is A. No reference FASTA is supplied to this workflow |
| Matching original 10x barcodes | Cell identifiers already stored in ARC and mgatk outputs; the workflow uses them to match mitochondrial measurements to retained cells. At least some retained cells must match; unmatched cells retain missing mgatk measurements |

Supply `--samples` with your sample CSV and `--stage` with the stage to analyse. For example, `--samples examples/samples.csv --stage E8.5` selects all rows labelled `E8.5`, including WT and mutant samples, and excludes other stages. Matching is case-sensitive; the workflow stops if no rows match.

Selected samples share a common peak set and are merged for stage-level RNA, ATAC and WNN analysis during the analysis phase. Outputs are written under `<outdir>/<stage>/`, such as `<outdir>/E8.5/`. Use the same `--stage` for review, cluster review and analysis.

Also supply the absolute MACS3 executable path with `--macs3_path` and, when using the Conda profile, the environment path with `--conda_env`. Cluster review and analysis additionally require your edited QC JSON through `--qc_json`. The commands below use the environment and MACS3 paths listed above.

## Run order: review, cluster review, analysis

Run all phases (review, cluster review, analysis) from the project directory with the same Nextflow work directory and stage. With `-resume`, unchanged upstream tasks can be reused. Use `-profile conda` for local execution, or `-profile slurm,conda` on Slurm. Pass site-specific queue/account settings with `-c site.config`.

First review raw QC distributions:

```bash
nextflow run main.nf \
  --samples /absolute/path/samples.csv --stage E8.5 \
  --outdir /absolute/path/stage_results \
  --conda_env /suffolk/WorkGenomicsD/yy605/conda_env/stage_multiome \
  --macs3_path /suffolk/WorkGenomicsD/yy605/conda_env/macs3/bin/macs3 \
  --phase review -profile slurm,conda -resume
```

Inspect `stage_results/E8.5/qc/review/`. Copy `qc_thresholds.template.json` to a separate file and enter your chosen numeric bounds. Its all-null thresholds deliberately fail validation: it is a review template, not a runnable configuration. Thresholds in `examples/qc_thresholds.example.json` are illustrative, not recommendations for your samples.

Then inspect clusters of the QC-passing cells before DoubletFinder:

```bash
nextflow run main.nf \
  --samples /absolute/path/samples.csv --stage E8.5 \
  --outdir /absolute/path/stage_results \
  --conda_env /suffolk/WorkGenomicsD/yy605/conda_env/stage_multiome \
  --macs3_path /suffolk/WorkGenomicsD/yy605/conda_env/macs3/bin/macs3 \
  --qc_json /absolute/path/edited_qc_thresholds.json \
  --phase cluster_review -profile slurm,conda -resume
```

Inspect `qc/cluster_review/`: per-cell cluster assignments, cluster-level QC summaries, violin plots and RNA markers. Add any low-quality cluster IDs to `doublet.exclude_clusters` for the corresponding sample. Empty exclusions mean no additional cluster removal; the workflow does not automatically identify low-quality biological populations. Complete cluster review before analysis: analysis now requires each sample's saved `<sample>.cluster_cells.tsv`. By default, it reads these files from `<outdir>/<stage>/qc/cluster_review/`; use `--cluster_review_dir /path/to/reviewed_tables` if they are stored elsewhere. If QC thresholds, clustering settings, input data or package versions change, review the resulting clusters again because their IDs may change.

Add cluster exclusions to a copy of the JSON used for cluster review, preserving its QC thresholds and clustering settings. Analysis uses the exact saved cluster labels, rather than recreating them; if you change inputs or QC/clustering settings, run cluster review again.

Finally run analysis with the reviewed JSON:

```bash
nextflow run main.nf \
  --samples /absolute/path/samples.csv --stage E8.5 \
  --outdir /absolute/path/stage_results \
  --conda_env /suffolk/WorkGenomicsD/yy605/conda_env/stage_multiome \
  --macs3_path /suffolk/WorkGenomicsD/yy605/conda_env/macs3/bin/macs3 \
  --qc_json /absolute/path/edited_qc_thresholds.json \
  --phase analysis -profile slurm,conda -resume
```

## QC and analysis JSON

The JSON has four sections: `thresholds`, optional `sample_overrides`, `clustering`, and optional `doublet`. An abridged settings example (retain all eight threshold entries from the template):

```json
{
  "clustering": {
    "algorithm": 3,
    "resolution": 2,
    "RNA_dims": "1:20",
    "RNA_umap_dims": "1:50",
    "ATAC_dims": "2:40",
    "ATAC_integrate_dims": "1:40",
    "integrate_rna": true,
    "integrate_atac": true,
    "min_cutoff": 10
  },
  "doublet": {
    "expected_rate": 0.05,
    "sample_rates": {},
    "exclude_clusters": {}
  }
}
```

For reviewed sample-specific decisions, use this structure; the cluster IDs and rate below are examples only:

```json
"doublet": {
  "expected_rate": 0.05,
  "sample_rates": {"5019_Hi_1": 0.06},
  "exclude_clusters": {"5019_Hi_1": ["3", "7"]}
}
```

Each metric must have `min` and `max`, with at least one numeric bound. `null` disables one side. Cells must pass all eight metrics; bounds are inclusive and non-finite values fail QC. `sample_overrides` replaces both bounds for each listed metric in a named sample. For example:

```json
"sample_overrides": {
  "5019_Hi_1": {"nCount_RNA": {"min": 1200, "max": 35000}}
}
```

Dimension strings are inclusive ranges. `RNA_dims` controls per-sample PCA/neighbours/clustering, DoubletFinder PCs, RNA CCA anchors, and RNA dimensions used by WNN. `RNA_umap_dims` controls both stage RNA UMAPs; stage PCA and integrated CCA retain enough components for this range. `ATAC_dims` controls ATAC neighbours, clustering, UMAP, reciprocal-LSI anchors and WNN. Its default drops LSI component 1. `ATAC_integrate_dims` controls corrected LSI components and must extend at least to the largest `ATAC_dims` component. `min_cutoff = 10` selects peaks with more than 10 total counts for LSI; it does not delete other peaks from the raw assay.

New templates enable both integrations. Each requires at least two samples. For a single-sample stage set both switches to `false`. Older JSON files without `integrate_rna` default to RNA integration; their explicit `integrate_atac` setting is retained. Missing `doublet` settings default to a 5% expected rate and no cluster exclusions. Choose the expected rate from your experiment's loading information; 5% is a pipeline default, not a universal biological estimate.

## Details for each step

Step labels follow the script filenames. Steps 03 and 04 are described together; 05a performs cluster review before 05b runs DoubletFinder during analysis. The Step 07–08 section covers modality combination and WNN; script 07 is a legacy helper, and the active workflow uses script 08.

### Step 01: call ATAC peaks per sample

Related tutorial: [Signac 1.2.0 multiome peak-calling tutorial](https://stuartlab.org/signac/1.2.0/articles/pbmc_multiomic#peak-calling-1).

`01_peaks.R` reads ARC filtered-cell peaks and fragments, builds a ChromatinAssay and calls MACS3 through Signac. A single grouping variable restricts peak calling to Cell Ranger called cells. It uses `GenomeInfoDb::keepStandardChromosomes(species = "Mus musculus", pruning.mode = "coarse")` to retain standard mouse chromosomes. It then removes entire peaks overlapping `Signac::blacklist_mm10` using `GenomicRanges::subsetByOverlaps(invert = TRUE)` before exporting the BED. When `chrM` is listed under `non_nuclear_contigs` in the reference used for the run, Cell Ranger ARC excludes mitochondrial fragments from the ATAC fragments file while retaining mitochondrial alignments in the ATAC BAM for analysis with mgatk. The effective mouse genome size is `1.87e9`.

Outputs: `STAGE/peaks/sample/<sample>.bed`. This uses all fragments belonging to the ARC filtered-cell barcodes, before the reviewed QC thresholds; it does not call peaks from non-cell barcodes or require a peak to be present in every sample.

### Step 02: create and recount a stage-level common peak set

`02_consensus_peaks.R` takes every sample peak BED within the chosen stage and uses `GenomicRanges::reduce` to create a union, merging overlapping or adjacent intervals. Default width filtering retains peaks with width strictly greater than 20 bp and strictly less than 10,000 bp (`--peak_min 20 --peak_max 10000`). These defaults are defined near the top of `main.nf` as `params.peak_min = 20` and `params.peak_max = 10000`, so you do not need to supply them on the command line. To override them, add arguments such as `--peak_min 50 --peak_max 5000` to your Nextflow command. Both boundaries are excluded: this example retains widths strictly greater than 50 bp and strictly less than 5,000 bp.

Output: `STAGE/peaks/common_peaks.bed`. The shared helper subsequently uses `FeatureMatrix` to recount these exact intervals in each sample's fragments, retaining the complete feature set even for absent chromosomes. Both review and filtering use these counts. This follows the [Signac merging procedure](https://stuartlab.org/signac/articles/merging): establish common intervals, quantify them separately, then merge.

### Step 03–04: review the eight RNA/ATAC QC metrics

`03_qc_metrics.R` creates RNA and common-peak ATAC objects for cells shared between modalities. It calculates the following metrics before applying user thresholds:

| Modality | Metric in files/JSON | Meaning |
|---|---|---|
| RNA | `nCount_RNA` | Total RNA UMI counts |
| RNA | `nFeature_RNA` | Number of detected RNA genes |
| RNA | `percent.mt` | RNA percentage from genes matching `^[Mm][Tt]-` |
| RNA | `percent.ribo` | RNA percentage from genes matching `^Rp[ls]` |
| ATAC | `nCount_ATAC` | Counts summed over the recounted common peaks |
| ATAC | `nFeature_ATAC` | Number of detected common peaks |
| ATAC | `nucleosome_signal` | Signac fragment-length nucleosome signal |
| ATAC | `TSS.enrichment` | Signac TSS enrichment using mm10 annotation |

Per-sample outputs are `<sample>.metrics.tsv` and `<sample>.distributions.pdf`, with one `ggpubr::ggviolin` page per metric and an inset boxplot. QC violin plots are generated by the shared `plot_qc_violins()` function in [scripts/00_common.R](scripts/00_common.R). Edit this function to change the theme or font sizes for per-sample, stage-level and cluster-review QC plots.

Both the per-sample `<sample>.metrics.tsv` files and the combined `all_cells_qc_metrics.tsv` are saved under `<outdir>/<stage>/qc/review/`, for example `results/E8.5/qc/review/`.

`04_qc_review.R` validates and combines the sample tables in sample-ID order, then writes:

| Output | Content |
|---|---|
| `all_cells_qc_metrics.tsv` | Combined per-cell measurements from all `<sample>.metrics.tsv` files in the selected stage; no additional per-cell metrics are calculated |
| `qc_summary.tsv` | One row per sample/metric, plus pooled `ALL_STAGE`: n, mean, SD, min, quantiles, median, max, scaled MAD and 3/5/7-MAD bounds |
| `qc_mad_candidates.tsv` | One row per sample/metric/multiplier with lower and upper bounds |
| `stage_qc_distributions.pdf` | Eight violin-plot pages comparing samples |
| `qc_thresholds.template.json` | Editable thresholds and analysis settings; written by `jsonlite::write_json()` in the task work directory, then copied to `<outdir>/<stage>/qc/review/qc_thresholds.template.json` |

Run `Rscript tests/qc_review.R` to check table combination, QC summaries, MAD bounds, template generation and plotting with synthetic data.

For k = 3, 5 or 7, bounds are `median ± k × mad(x, constant = 1.4826)` and appear as `mad_k_lower`/`mad_k_upper` in the summary. Non-finite values are excluded from statistics. Bounds are not clipped to zero and are not automatically applied. Zero MAD gives equal lower/upper bounds; inspect the distribution before choosing thresholds.

### Step 05a: review low-quality clusters before DoubletFinder

`05a_review_clusters.R`, run by `--phase cluster_review`, applies the reviewed eight-metric QC thresholds to each sample. It log-normalizes RNA, selects 2,000 variable genes, scales, computes individual PCA, finds neighbours with `RNA_dims`, and clusters with the configured algorithm/resolution. The saved `<sample>.cluster_cells.tsv` contains the barcode-to-cluster assignments that step `05b` uses for cluster removal.

Outputs in `STAGE/qc/cluster_review/`: `<sample>.cluster_cells.tsv`, `<sample>.cluster_summary.tsv`, `<sample>.cluster_qc.pdf`, and `<sample>.cluster_markers.tsv`. Marker tables retain only rows with finite `p_val_adj < 0.05`. Review RNA counts, mitochondrial percentages and marker informativeness together. Enter unwanted cluster IDs as strings in `doublet.exclude_clusters`. This provides the manual review/removal route described by the [DoubletFinder authors](https://github.com/chris-mcginnis-ucsf/DoubletFinder#best-practices-for-scrna-seq-data-generated-without-sample-multiplexing). For a related example of reviewing QC metrics by cluster, removing low-quality clusters and repeating preprocessing, see the [snRNA-seq Guide: cluster QC](https://lahuuki.github.io/rstats/2026-03-19-snRNAseqGuide/index.html).

### Step 05b: joint QC, cluster removal and DoubletFinder

`05b_qc_doubletfinder.R` applies all eight QC thresholds and reads the saved `<sample>.cluster_cells.tsv` from step 05a. It validates the sample ID, unique barcodes, cluster labels, QC-passing cell set and QC measurements, aligns labels by barcode, and removes the listed clusters without repeating the preliminary clustering. It then computes normalization, variable features, PCA, neighbours and clusters on the retained cells before running DoubletFinder. Unknown exclusion IDs, mismatched cells and changed QC measurements fail explicitly; rerun cluster review when inputs or QC decisions change. The original reviewed labels remain in the RNA metadata as `reviewed_qc_cluster`; newly computed `pre_doublet_clusters` are used for DoubletFinder's homotypic correction. At least 100 cells must remain before DoubletFinder.

DoubletFinder runs separately on each sample, using log-normalized RNA (`sct = FALSE`) and the configured RNA dimensions. A pN/pK sweep selects the pK with maximum BCmetric. The expected doublet count uses the configured sample rate and a homotypic correction derived from the reclustered RNA populations. Inspect the estimates for your experiment; this correction uses clusters as a proxy for cell types. At least 50 singlets must remain.

The resulting singlet barcodes select both RNA and ATAC. Cell names gain the sample prefix; metadata retains `sample_id` and the original `barcode` for mgatk matching. Outputs are individual `<sample>.rna.qc.rds` and `<sample>.atac.qc.rds` checkpoints in the Nextflow work directory, and a published `<sample>.qc.tsv` recording raw/QC/cluster/singlet counts, excluded cluster IDs, pK, rate and expected doublets.

Cell counts are reported per sample in `<outdir>/<stage>/qc/filtered/<sample>.qc.tsv` and in the task's `.command.err` log. `raw_joint` is the number of starting cells with both modalities. Each removal count is relative to the immediately preceding step; if no clusters are excluded, `cluster_removed` is zero.

| Step | Cells remaining | Cells removed |
| --- | --- | --- |
| 1. Initial QC metric filtering | `qc_pass` | `qc_removed` |
| 2. Reviewed-cluster removal, before reclustering | `cluster_pass` | `cluster_removed` |
| 3. DoubletFinder singlet selection | `df_singlets` | `df_removed` |

### Step 06a: individual RNA preprocessing, then stage merge and integration

`06a_preprocess_rna.R` receives the JSON and all filtered RNA RDS files for the stage. For each sample, in memory, it performs log-normalization, 2,000-variable-gene selection, scaling, individual PCA, neighbours and sample clustering. It then adds:

- `CellCycleScoring`: `S.Score`, `G2M.Score` and `Phase`, using case-insensitive matches of Seurat's 2019 marker panel to observed mouse symbols. This is symbol matching, not the online gprofiler orthology conversion used in YY_04b. At least five markers per panel must match. No cell-cycle regression is applied.
- `cellXY::classifySex(counts, genome = "Mm", qc = FALSE)`: stored as `sex_prediction`. The cells already passed pipeline QC; uncertain predictions remain missing. These are expression-based predictions, not sample sex assignments. cellXY internally converts counts to a dense matrix, so large samples may require substantial RAM.
- SCT v2 normalization with 2,000 variable features per sample.

After individual processing it selects shared SCT integration features, prepares SCT residuals, and merges samples. It retains joined log-normalized RNA for marker testing, computes `rna_merged_pca` on SCT features, and creates `rna_unintegrated_clusters` and `rna_unintegrated_umap`. When enabled, SCT CCA `IntegrateLayers` creates `integrated.cca`, followed by `rna_integrated_clusters`, `rna_integrated_umap`, and `PrepSCTFindMarkers`. This implements the YY_05 merge/integration sequence with explicit dimension settings; see [Seurat CCA integration](https://satijalab.org/seurat/reference/ccaintegration).

Published outputs: `rna_unintegrated.rds`, optional `rna_integrated.rds`, selected `rna_stage.rds`, and `rna_cell_annotations.tsv` containing per-cell sample, cycle scores/phase and sex prediction. No individual preprocessed RDS files are saved.

### Step 06b: individual ATAC preprocessing, then stage merge and integration

`06b_preprocess_atac.R` receives the JSON and all filtered ATAC RDS files. In memory it runs sample-level `FindTopFeatures`, TF-IDF and SVD/LSI; these sample LSI reductions are required for reciprocal-LSI anchors. It merges the objects while retaining fragment mappings, then recomputes stage TF-IDF/LSI and generates `atac_unintegrated_clusters` and `atac_unintegrated_umap`.

When enabled, `FindIntegrationAnchors(reduction = "rlsi")` and `IntegrateEmbeddings` create `integrated_lsi`, followed by `atac_integrated_clusters` and `atac_integrated_umap`. Both integration steps use a weighting neighbourhood capped by the smallest sample size. This follows YY_06 and the [Signac ATAC integration workflow](https://stuartlab.org/signac/articles/integrate_atac.html).

Published outputs: `atac_unintegrated.rds`, optional `atac_integrated.rds`, and selected `atac_stage.rds`. No individual preprocessed RDS files are saved. Compare integrated and unintegrated results because integration can reduce sample-associated biological variation.

### Step 06c: read precomputed mgatk results per sample

`06c_mgatk_metrics.R` reads completed mgatk results with `ReadMGATK`; it does not run mitochondrial alignment, variant calling or mgatk itself. It matches the original retained RNA barcodes, checks reference A at m.5019, and calculates total and strand-specific A/G counts, locus coverage, mtDNA sequencing depth and G/(A+G) heteroplasmy. Missing mgatk barcodes are retained with missing measurements.

Output: `STAGE/mgatk/<sample>.mtDNA.tsv`. `hf_5019_covered` is missing below 10 A+G reads; `depth_5019` includes all four bases. `mtDNA_depth` measures sequencing coverage, not absolute mtDNA copy number.

### Step 07–08: combine modalities and build WNN

`08_wnn_stage.R` aligns the stage RNA and ATAC cells, verifies their sample assignments, and attaches ATAC assays, QC metadata, reductions and modality cluster labels to the RNA object. It uses `integrated.cca` and `integrated_lsi` when their switches are enabled, or `rna_merged_pca` and `lsi` otherwise. It preserves stage modality UMAPs and cycle/sex metadata.

`FindMultiModalNeighbors` uses `RNA_dims` and `ATAC_dims`, storing modality weights. WNN UMAP uses the weighted neighbour graph; clusters use the `wsnn` graph and configured algorithm/resolution. Selected reductions are recorded in `obj@misc$wnn_reductions`.

Output: `STAGE/objects/wnn.rds`, including `wnn.umap` and `wsnn_clusters`. The old `07_combine_modalities.R` remains a legacy helper and is not invoked by Nextflow.

### Step 09: attach mgatk measurements to the WNN object

`09_add_mgatk.R` combines the sample metric tables, rejects duplicated/mismatched cell IDs, orders measurements to the WNN cell order, and attaches all measured mtDNA fields, including strand-specific heteroplasmy.

Outputs: `STAGE/objects/wnn_mgatk.rds` and `all_cells_mtDNA.tsv`.

### Step 10: RNA markers for each WNN cluster

`10_cluster_markers.R` joins RNA layers, uses `wsnn_clusters` as identities, and calls `FindAllMarkers` on the RNA assay. It compares each cluster with all other cells using Wilcoxon tests, `only.pos = TRUE`, `min.pct = 0.1` and `logfc.threshold = 0.25`. The output retains only rows with finite `p_val_adj < 0.05`.

Output: `STAGE/deg/cluster_markers.tsv`. These are positive RNA cluster markers; this step does not perform ATAC differential accessibility, sample-level pseudobulk tests or Hi/Lo comparisons. Cell-level p-values do not establish embryo-level differential expression.

## Optional: estimate mtDNA copy number from mgatk and ARC

`scripts/11_estimate_mtdna_copy_number.R` is a standalone, commented script; it is not invoked by `main.nf`. It estimates mtDNA copies as `ploidy × mitochondrial mean depth / nuclear mean depth`, with a default nuclear ploidy of 2. It reads genome-wide mitochondrial depth from mgatk's `*.depthTable.txt` and measures nuclear sequenced-base coverage from the ARC ATAC BAM. It does not use heteroplasmy or coverage at m.5019 as a copy-number measurement. See the [mgatk output definitions](https://github.com/caleblareau/mgatk/wiki/Outputs).

First export original barcodes for one sample from its published mgatk metrics table:

```r
metrics <- read.delim("/path/to/stage_results/E8.5/mgatk/5019_Hi_1.mtDNA.tsv")
writeLines(as.character(metrics$barcode), "5019_Hi_1.barcodes.txt")
```

Then run in the existing R environment:

```bash
conda run -p /suffolk/WorkGenomicsD/yy605/conda_env/stage_multiome \
  Rscript scripts/11_estimate_mtdna_copy_number.R \
  --sample 5019_Hi_1 \
  --bam /path/to/arc/outs/atac_possorted_bam.bam \
  --fragments /path/to/arc/outs/atac_fragments.tsv.gz \
  --barcode-map /path/to/arc/outs/per_barcode_metrics.csv \
  --barcodes 5019_Hi_1.barcodes.txt \
  --mgatk-dir /path/to/mgatk/final \
  --blacklist /path/to/mm10-blacklist.bed \
  --mapq 30 --baseq 20 --keep-duplicates false \
  --out-prefix /path/to/copy_number/5019_Hi_1
```

`--barcode-map` is optional when BAM CB tags and mgatk identifiers already match the original ARC barcodes. Otherwise it translates `atac_barcode` to the canonical `barcode` using ARC's metrics table, for both BAM and mgatk matching. Inputs are for one sample per invocation. The BAM must have a valid reference header; chromosome lengths are read from it. Use an unaggregated BAM from the same ATAC library used for mgatk. Required packages are Rsamtools, GenomicRanges, GenomicAlignments and their dependencies, already listed in `environment.yml`.

The script streams fragments, scores their endpoints in 500 bp bins on chr1–chr19, excludes bins overlapping the supplied mm10 blacklist, and selects scores strictly above the 98.5th percentile. It excludes sex chromosomes and mitochondria from the nuclear reference. It ignores fragment readSupport (PCR duplicate multiplicity). Reference selection can be adjusted with `--bin-size` and `--quantile`; zero-score eligible bins participate in the percentile. It uses no neighbouring-bin expansion. BAM coverage includes only A/C/G/T bases passing quality filters inside the exact selected regions, respecting insertions, deletions and soft clipping. Both mates contribute separately if they overlap. The BAM is scanned in chunks; per-read R counting can take substantial time on a full library, so start with a small pilot.

For comparisons across samples, use a fixed, representative nuclear reference set. Supply it with `--reference-bed /path/to/shared_reference.bed` instead of `--fragments` and `--blacklist`. In this mode the BED is used as supplied; it must already exclude blacklisted regions and contain only mouse autosomes. The script reduces overlapping intervals so reference length and coverage are not counted twice. A reference selected from one sample may favour that sample's cell composition; assess accessibility stability across the cell types and samples being compared.

Outputs for each prefix:

| File | Content |
|---|---|
| `.mtDNA_copy_number.tsv` | Original barcodes, nuclear aligned bases/read counts, nuclear reference length, both mean depths, assumed ploidy, raw estimate, estimate passing coverage checks, and QC status |
| `.nuclear_reference.bed` | Exact selected nuclear intervals, reusable with `--reference-bed` |
| `.settings.txt` | Arguments, reference length, mgatk depth-table path, mate-overlap policy and R/package versions |

Defaults flag cells below 1,000 nuclear aligned bases or mitochondrial mean depth 10; these are configurable coverage checks, not validated universal cutoffs. Missing mgatk measurements and zero nuclear coverage remain missing estimates. Low-coverage cells retain their raw estimate for inspection, but `mtDNA_copy_estimate` is missing unless `qc_status` is `pass`. Use `--min-nuclear-bases` and `--min-mt-depth` to adjust these checks.

The script uses the supplied genome-wide mgatk mean depth without applying the notebook's human mitochondrial length or deletion correction. Verify mitochondrial reference and coverage exclusions in your mgatk run. **Match the nuclear mapping/base-quality thresholds, duplicate handling and mate-overlap conventions to the actual mgatk run.** The script cannot recover those settings from the depth table. Different filters or libraries can distort the coverage ratio. The diploid, equal-sampling assumption can also vary with cell cycle, nuclear copy-number changes, accessibility and mitochondrial recovery; interpret the result as an ATAC-derived estimate and validate absolute copy numbers independently when needed.

Run `Rscript tests/mtdna_copy_number.R` for synthetic CIGAR, base-quality, duplicate-filtering, barcode-translation, missing-data, reference-selection and BAM-to-estimate checks. These checks do not validate the biological calibration.

## Workflow files and execution

`main.nf` wires every active script and stages files between tasks. `nextflow.config` supplies local/Conda/Slurm execution settings. `scripts/00_common.R` contains shared import, QC, validation, clustering, scoring and annotation helpers; it is not a standalone process. `examples/` contains the manifest and illustrative JSON.

RNA and ATAC stage processing can run concurrently, while mgatk extraction runs per sample. WNN waits for both stage objects; mgatk attachment waits for WNN and all mgatk tables; marker detection waits for the annotated final object. The stage RNA/ATAC tasks each request 96 GB; adjust process resources in a site config if needed.

Run the dependency-light regression checks with `Rscript tests/regression.R`. These parse every R script and exercise configuration validation, inclusive QC boundaries, non-finite QC rejection, and missing mgatk barcode handling using synthetic data. They do not exercise the full scientific analysis. This workspace has Java 8 and lacks Signac, DoubletFinder and cellXY, so full workflow execution requires the production environment and a small real stage.

Samples are validated before tasks launch and processed in sample-ID order. Shell arguments support paths containing spaces and apostrophes. An R-source checksum is included in task commands, so editing any R script invalidates cached tasks on `-resume`. ARC and mgatk directories remain shared absolute paths: changing their contents in place is not detected by Nextflow; use a fresh work directory when replacing those data. Package/environment changes also require a fresh work directory.
