# Per-stage m.5019A>G Multiome workflow

This Nextflow workflow analyses one mouse developmental stage at a time (mm10). It calls peaks from each sample's ATAC fragments, creates one common peak set for the stage, reviews eight RNA/ATAC QC metrics, removes reviewed low-quality clusters and RNA doublets per sample, preprocesses and integrates RNA and ATAC, builds WNN clusters, attaches precomputed mgatk results, and identifies RNA markers for each WNN cluster.

RNA and ATAC always use the same retained singlet barcodes. Per-sample preprocessing takes place in memory inside the combined `06a` and `06b` stage tasks; these tasks save stage objects only. Step `05` still saves individual QC/singlet objects as workflow inputs and checkpoints.

## Requirements and inputs

Use Nextflow with Java 17 or later and a shared R environment containing Seurat 5, Signac 1.x, Matrix, sctransform, DoubletFinder (Seurat v5 compatible), GenomicRanges, GenomeInfoDb, rtracklayer, jsonlite, ggpubr, ggplot2, EnsDb.Mmusculus.v79 and cellXY. `glmGamPoi` is recommended for SCT fitting; without it, the script uses the slower Poisson implementation. SLM clustering uses algorithm 3; Leiden (algorithm 4) requires its additional dependencies.

Install MACS3 separately and supply its absolute executable path. Example installation in an activated Python environment:

```bash
python -m pip install macs3
macs3 --version
```

Check the [MACS3 installation documentation](https://macs3-project.github.io/MACS/docs/INSTALL.html) for the Python version and compiler requirements of the release you install. The workflow does not install software automatically. Signac 1.x accepts the MACS3 executable through its `macs2.path` argument.

Example R package installation, adapted to your existing environment:

```r
install.packages(c("Seurat", "Signac", "jsonlite", "ggpubr", "remotes"))
BiocManager::install(c("GenomicRanges", "GenomeInfoDb", "rtracklayer",
                      "EnsDb.Mmusculus.v79", "glmGamPoi"))
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
remotes::install_github("phipsonlab/cellXY")
```

Use [cellXY](https://github.com/phipsonlab/cellXY) for `classifySex`; current Bioconductor speckle releases do not provide the older sex-classification interface. Record/pin package versions in the environment used for production runs.

The sample CSV must contain `stage,sample_id,arc_out,mgatk_dir`; see `examples/samples.csv`. This project currently accepts sample IDs starting with `5019_` and letters, numbers, `_`, `.` or `-`. Each ID must be unique within its stage. All execution nodes must be able to read the absolute input paths:

| Input | Required content |
|---|---|
| `arc_out` | Cell Ranger ARC `outs` directory containing `filtered_feature_bc_matrix.h5`, `atac_fragments.tsv.gz` and its `.tbi` index |
| RNA features | Mouse gene symbols; cell-cycle scoring and sex inference require symbols rather than Ensembl IDs |
| `mgatk_dir` | Completed mgatk output readable by `Signac::ReadMGATK`, usually the `final` directory |
| Reference | mm10 for ARC and annotation; mitochondrial position 5019 must have reference allele A |
| Barcodes | Original 10x barcodes must agree between ARC and mgatk |

## Run order: review, cluster review, analysis

Run all phases from the project directory with the same Nextflow work directory and stage. With `-resume`, unchanged upstream tasks can be reused. Use `-profile conda` for local execution, or `-profile slurm,conda` on Slurm. Pass site-specific queue/account settings with `-c site.config`.

First review raw QC distributions:

```bash
nextflow run main.nf \
  --samples /absolute/path/samples.csv --stage E8.5 \
  --outdir /absolute/path/stage_results \
  --conda_env /shared/path/to/r-environment \
  --macs3_path /shared/path/to/macs3 \
  --phase review -profile slurm,conda -resume
```

Inspect `stage_results/E8.5/qc/review/`. Copy `qc_thresholds.template.json` to a separate file and enter your chosen numeric bounds. Its all-null thresholds deliberately fail validation: it is a review template, not a runnable configuration. Thresholds in `examples/qc_thresholds.example.json` are illustrative, not recommendations for your samples.

Then inspect clusters of the QC-passing cells before DoubletFinder:

```bash
nextflow run main.nf \
  --samples /absolute/path/samples.csv --stage E8.5 \
  --outdir /absolute/path/stage_results \
  --conda_env /shared/path/to/r-environment \
  --macs3_path /shared/path/to/macs3 \
  --qc_json /absolute/path/edited_qc_thresholds.json \
  --phase cluster_review -profile slurm,conda -resume
```

Inspect `qc/cluster_review/`: per-cell cluster assignments, cluster-level QC summaries, violin plots and RNA markers. Add any low-quality cluster IDs to `doublet.exclude_clusters` for the corresponding sample. Empty exclusions mean no additional cluster removal; the workflow does not automatically identify low-quality biological populations. Cluster review is a separate phase and is not automatically enforced by `analysis`. To follow the full checklist, complete it before analysis. If QC thresholds, clustering settings, input data or package versions change, review the resulting clusters again because their IDs may change.

Finally run analysis with the reviewed JSON:

```bash
nextflow run main.nf \
  --samples /absolute/path/samples.csv --stage E8.5 \
  --outdir /absolute/path/stage_results \
  --conda_env /shared/path/to/r-environment \
  --macs3_path /shared/path/to/macs3 \
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

### 01: call ATAC peaks per sample

`01_peaks.R` reads ARC filtered-cell peaks and fragments, builds a ChromatinAssay and calls MACS3 through Signac. A single grouping variable restricts peak calling to Cell Ranger called cells. It retains nuclear `chr1`–`chr19`, `chrX` and `chrY`-style chromosome names (the implementation accepts numeric chromosome names generally). The effective mouse genome size is `1.87e9`.

Outputs: `STAGE/peaks/sample/<sample>.bed`. This uses all fragments belonging to the ARC filtered-cell barcodes, before the reviewed QC thresholds; it does not call peaks from non-cell barcodes or require a peak to be present in every sample.

### 02: create and recount a stage-level common peak set

`02_consensus_peaks.R` takes every sample peak BED within the chosen stage and uses `GenomicRanges::reduce` to create a union, merging overlapping or adjacent intervals. Default width filtering retains peaks with width strictly greater than 20 bp and strictly less than 10,000 bp (`--peak_min 20 --peak_max 10000`).

Output: `STAGE/peaks/common_peaks.bed`. The shared helper subsequently uses `FeatureMatrix` to recount these exact intervals in each sample's fragments, retaining the complete feature set even for absent chromosomes. Both review and filtering use these counts. This follows the [Signac merging procedure](https://stuartlab.org/signac/articles/merging): establish common intervals, quantify them separately, then merge. The ATAC stage step rejects inputs with different peak names/order.

### 03–04: review the eight RNA/ATAC QC metrics

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

Per-sample outputs are `<sample>.metrics.tsv` and `<sample>.distributions.pdf`, with one `ggpubr::ggviolin` page per metric and an inset boxplot. No MAD lines are added to these plots. ATAC nCount is the common-peak count, not ARC's total fragment count.

`04_qc_review.R` collects the sample tables and writes:

| Output | Content |
|---|---|
| `all_cells_qc_metrics.tsv` | All samples' per-cell metrics |
| `qc_summary.tsv` | One row per sample/metric, plus pooled `ALL_STAGE`: n, mean, SD, min, quantiles, median, max, scaled MAD and 3/5/7-MAD bounds |
| `qc_mad_candidates.tsv` | One row per sample/metric/multiplier with lower and upper bounds |
| `stage_qc_distributions.pdf` | Eight violin-plot pages comparing samples |
| `qc_thresholds.template.json` | Editable thresholds and analysis settings |

For k = 3, 5 or 7, bounds are `median ± k × mad(x, constant = 1.4826)` and appear as `mad_k_lower`/`mad_k_upper` in the summary. Non-finite values are excluded from statistics. Bounds are not clipped to zero and are not automatically applied. Zero MAD gives equal lower/upper bounds; inspect the distribution before choosing thresholds.

### 05a: review low-quality clusters before DoubletFinder

`05a_review_clusters.R`, run by `--phase cluster_review`, applies the reviewed eight-metric QC thresholds to each sample. It log-normalizes RNA, selects 2,000 variable genes, scales, computes individual PCA, finds neighbours with `RNA_dims`, and clusters with the configured algorithm/resolution. The same deterministic helper is used by step `05` to reproduce the pre-doublet cluster labels.

Outputs in `STAGE/qc/cluster_review/`: `<sample>.cluster_cells.tsv`, `<sample>.cluster_summary.tsv`, `<sample>.cluster_qc.pdf`, and `<sample>.cluster_markers.tsv`. Review RNA counts, mitochondrial percentages and marker informativeness together. Enter unwanted cluster IDs as strings in `doublet.exclude_clusters`. This provides the manual review/removal route described by the [DoubletFinder authors](https://github.com/chris-mcginnis-ucsf/DoubletFinder#best-practices-for-scrna-seq-data-generated-without-sample-multiplexing).

### 05: joint QC, cluster removal and DoubletFinder

`05_qc_doubletfinder.R` applies all eight QC thresholds, recreates the pre-doublet RNA clusters, removes the listed clusters, and recomputes normalization, variable features, PCA, neighbours and clusters before running DoubletFinder. Unknown exclusion IDs fail explicitly. At least 100 cells must remain before DoubletFinder.

DoubletFinder runs separately on each sample, using log-normalized RNA (`sct = FALSE`) and the configured RNA dimensions. A pN/pK sweep selects the pK with maximum BCmetric. The expected doublet count uses the configured sample rate and a homotypic correction derived from the reclustered RNA populations. Inspect the estimates for your experiment; this correction uses clusters as a proxy for cell types. At least 50 singlets must remain.

The resulting singlet barcodes select both RNA and ATAC. Cell names gain the sample prefix; metadata retains `sample_id` and the original `barcode` for mgatk matching. Outputs are individual `<sample>.rna.qc.rds` and `<sample>.atac.qc.rds` checkpoints in the Nextflow work directory, and a published `<sample>.qc.tsv` recording raw/QC/cluster/singlet counts, excluded cluster IDs, pK, rate and expected doublets.

### 06a: individual RNA preprocessing, then stage merge and integration

`06a_preprocess_rna.R` receives the JSON and all filtered RNA RDS files for the stage. For each sample, in memory, it performs log-normalization, 2,000-variable-gene selection, scaling, individual PCA, neighbours and sample clustering. It then adds:

- `CellCycleScoring`: `S.Score`, `G2M.Score` and `Phase`, using case-insensitive matches of Seurat's 2019 marker panel to observed mouse symbols. This is symbol matching, not the online gprofiler orthology conversion used in YY_04b. At least five markers per panel must match. No cell-cycle regression is applied.
- `cellXY::classifySex(counts, genome = "Mm", qc = FALSE)`: stored as `sex_prediction`. The cells already passed pipeline QC; uncertain predictions remain missing. These are expression-based predictions, not sample sex assignments. cellXY internally converts counts to a dense matrix, so large samples may require substantial RAM.
- SCT v2 normalization with 2,000 variable features per sample.

After individual processing it selects shared SCT integration features, prepares SCT residuals, and merges samples. It retains joined log-normalized RNA for marker testing, computes `rna_merged_pca` on SCT features, and creates `rna_unintegrated_clusters` and `rna_unintegrated_umap`. When enabled, SCT CCA `IntegrateLayers` creates `integrated.cca`, followed by `rna_integrated_clusters`, `rna_integrated_umap`, and `PrepSCTFindMarkers`. This implements the YY_05 merge/integration sequence with explicit dimension settings; see [Seurat CCA integration](https://satijalab.org/seurat/reference/ccaintegration).

Published outputs: `rna_unintegrated.rds`, optional `rna_integrated.rds`, selected `rna_stage.rds`, and `rna_cell_annotations.tsv` containing per-cell sample, cycle scores/phase and sex prediction. No individual preprocessed RDS files are saved.

### 06b: individual ATAC preprocessing, then stage merge and integration

`06b_preprocess_atac.R` receives the JSON and all filtered ATAC RDS files. In memory it runs sample-level `FindTopFeatures`, TF-IDF and SVD/LSI; these sample LSI reductions are required for reciprocal-LSI anchors. It merges the objects while retaining fragment mappings, then recomputes stage TF-IDF/LSI and generates `atac_unintegrated_clusters` and `atac_unintegrated_umap`.

When enabled, `FindIntegrationAnchors(reduction = "rlsi")` and `IntegrateEmbeddings` create `integrated_lsi`, followed by `atac_integrated_clusters` and `atac_integrated_umap`. Both integration steps use a weighting neighbourhood capped by the smallest sample size. This follows YY_06 and the [Signac ATAC integration workflow](https://stuartlab.org/signac/articles/integrate_atac.html).

Published outputs: `atac_unintegrated.rds`, optional `atac_integrated.rds`, and selected `atac_stage.rds`. No individual preprocessed RDS files are saved. Compare integrated and unintegrated results because integration can reduce sample-associated biological variation.

### 06c: read precomputed mgatk results per sample

`06c_mgatk_metrics.R` reads completed mgatk results with `ReadMGATK`; it does not run mitochondrial alignment, variant calling or mgatk itself. It matches the original retained RNA barcodes, checks reference A at m.5019, and calculates total and strand-specific A/G counts, locus coverage, mtDNA sequencing depth and G/(A+G) heteroplasmy. Missing mgatk barcodes are retained with missing measurements.

Output: `STAGE/mgatk/<sample>.mtDNA.tsv`. `hf_5019_covered` is missing below 10 A+G reads; `depth_5019` includes all four bases. `mtDNA_depth` measures sequencing coverage, not absolute mtDNA copy number.

### 08: combine modalities and build WNN

`08_wnn_stage.R` aligns the stage RNA and ATAC cells, verifies their sample assignments, and attaches ATAC assays, QC metadata, reductions and modality cluster labels to the RNA object. It uses `integrated.cca` and `integrated_lsi` when their switches are enabled, or `rna_merged_pca` and `lsi` otherwise. It preserves stage modality UMAPs and cycle/sex metadata.

`FindMultiModalNeighbors` uses `RNA_dims` and `ATAC_dims`, storing modality weights. WNN UMAP uses the weighted neighbour graph; clusters use the `wsnn` graph and configured algorithm/resolution. Selected reductions are recorded in `obj@misc$wnn_reductions`.

Output: `STAGE/objects/wnn.rds`, including `wnn.umap` and `wsnn_clusters`. The old `07_combine_modalities.R` remains a legacy helper and is not invoked by Nextflow.

### 09: attach mgatk measurements to the WNN object

`09_add_mgatk.R` combines the sample metric tables, rejects duplicated/mismatched cell IDs, orders measurements to the WNN cell order, and attaches all measured mtDNA fields, including strand-specific heteroplasmy.

Outputs: `STAGE/objects/wnn_mgatk.rds` and `all_cells_mtDNA.tsv`.

### 10: RNA markers for each WNN cluster

`10_cluster_markers.R` joins RNA layers, uses `wsnn_clusters` as identities, and calls `FindAllMarkers` on the RNA assay. It compares each cluster with all other cells using Wilcoxon tests, `only.pos = TRUE`, `min.pct = 0.1` and `logfc.threshold = 0.25`.

Output: `STAGE/deg/cluster_markers.tsv`. These are positive RNA cluster markers; this step does not perform ATAC differential accessibility, sample-level pseudobulk tests or Hi/Lo comparisons. Cell-level p-values do not establish embryo-level differential expression.

## Workflow files and execution

`main.nf` wires every active script and stages files between tasks. `nextflow.config` supplies local/Conda/Slurm execution settings. `scripts/00_common.R` contains shared import, QC, validation, clustering, scoring and annotation helpers; it is not a standalone process. `examples/` contains the manifest and illustrative JSON.

RNA and ATAC stage processing can run concurrently, while mgatk extraction runs per sample. WNN waits for both stage objects; mgatk attachment waits for WNN and all mgatk tables; marker detection waits for the annotated final object. The stage RNA/ATAC tasks each request 96 GB; adjust process resources in a site config if needed.

The source audit covers all scripts, workflow wiring, example files and this README. Synthetic QC, RNA integration and WNN tests have been run in this workspace. Full cluster-review/DoubletFinder, sex-inference, ATAC integration and cluster execution require the missing local dependencies (Signac, DoubletFinder, cellXY and Java 17+); validate a small real stage in the production environment before a full run.
