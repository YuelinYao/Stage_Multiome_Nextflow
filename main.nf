nextflow.enable.dsl = 2

params.samples = null
params.stage = null
params.outdir = 'results'
params.macs3_path = null
params.qc_json = null
params.phase = 'review'
params.peak_min = 20
params.peak_max = 10000

if (!params.samples || !params.stage) error 'Provide --samples manifest.csv and --stage E8.5'
if (!(params.phase in ['review', 'cluster_review', 'analysis'])) error '--phase must be review, cluster_review or analysis'
if (!params.macs3_path) error 'Provide --macs3_path /absolute/path/to/macs3'
if (!file(params.macs3_path).isAbsolute()) error '--macs3_path must be absolute'
if (!java.nio.file.Files.isRegularFile(file(params.macs3_path)) ||
    !java.nio.file.Files.isExecutable(file(params.macs3_path)))
  error "MACS3 executable not found or not executable: ${params.macs3_path}"
if (params.phase in ['cluster_review', 'analysis'] && (!params.qc_json || !file(params.qc_json).exists()))
  error 'Cluster review and analysis require --qc_json /absolute/path/to/edited_qc_thresholds.json'
if (params.enable_conda) {
  if (!params.conda_env) error 'The conda profile requires --conda_env /shared/path/to/env (or environment.yml)'
  if (!file(params.conda_env).exists()) error "Conda environment or YAML not found: ${params.conda_env}"
}

process PEAK_SAMPLE {
  tag "$id"
  cpus 4
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/peaks/sample", mode: 'copy'
  input:
  tuple val(id), val(arc), val(mgatk)
  output:
  tuple val(id), path("${id}.bed")
  script:
  """
  Rscript ${projectDir}/scripts/01_peaks.R '${id}' '${arc}' '${params.macs3_path}' '${id}.bed'
  """
}

process COMMON_PEAKS {
  cpus 2
  memory '12 GB'
  publishDir "${params.outdir}/${params.stage}/peaks", mode: 'copy'
  input:
  path beds
  output:
  path 'common_peaks.bed'
  script:
  """
  Rscript ${projectDir}/scripts/02_consensus_peaks.R common_peaks.bed ${params.peak_min} ${params.peak_max} ${beds.join(' ')}
  """
}

process QC_METRICS {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/review", mode: 'copy'
  input:
  tuple val(id), val(arc), val(mgatk), path(peaks)
  output:
  tuple val(id), path("${id}.metrics.tsv"), path("${id}.distributions.pdf")
  script:
  """
  Rscript ${projectDir}/scripts/03_qc_metrics.R '${id}' '${arc}' '${peaks}' '${id}.metrics.tsv' '${id}.distributions.pdf'
  """
}

process QC_REVIEW {
  cpus 2
  memory '12 GB'
  publishDir "${params.outdir}/${params.stage}/qc/review", mode: 'copy'
  input:
  path metric_files
  output:
  path 'qc_summary.tsv'
  path 'qc_mad_candidates.tsv'
  path 'qc_thresholds.template.json'
  path 'all_cells_qc_metrics.tsv'
  path 'stage_qc_distributions.pdf'
  script:
  """
  Rscript ${projectDir}/scripts/04_qc_review.R qc_summary.tsv qc_mad_candidates.tsv qc_thresholds.template.json stage_qc_distributions.pdf ${metric_files.join(' ')}
  """
}

process PRE_DOUBLET_REVIEW {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/cluster_review", mode: 'copy'
  input:
  tuple val(id), val(arc), val(mgatk), path(peaks), path(json)
  output:
  path "${id}.cluster_cells.tsv"
  path "${id}.cluster_summary.tsv"
  path "${id}.cluster_qc.pdf"
  path "${id}.cluster_markers.tsv"
  script:
  """
  Rscript ${projectDir}/scripts/05a_review_clusters.R '${id}' '${arc}' '${peaks}' '${json}'
  """
}

process QC_FILTER_DOUBLET {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/filtered", mode: 'copy', pattern: '*.qc.tsv'
  input:
  tuple val(id), val(arc), val(mgatk), path(peaks), path(json)
  output:
  tuple val(id), val(mgatk), path("${id}.rna.qc.rds"), path("${id}.atac.qc.rds"), path("${id}.qc.tsv")
  script:
  """
  Rscript ${projectDir}/scripts/05_qc_doubletfinder.R '${id}' '${arc}' '${peaks}' '${json}' '${id}.rna.qc.rds' '${id}.atac.qc.rds' '${id}.qc.tsv'
  """
}

process MERGE_INTEGRATE_RNA {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: 'rna_*'
  input:
  path objects
  path json
  output:
  path 'rna_stage.rds', emit: stage
  path 'rna_cell_annotations.tsv', emit: annotations
  path 'rna_unintegrated.rds', emit: unintegrated
  path 'rna_integrated.rds', optional: true, emit: integrated
  script:
  """
  Rscript ${projectDir}/scripts/06a_preprocess_rna.R '${json}' ${objects.join(' ')}
  """
}

process MERGE_INTEGRATE_ATAC {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: 'atac_*.rds'
  input:
  path objects
  path json
  output:
  path 'atac_stage.rds', emit: stage
  path 'atac_unintegrated.rds', emit: unintegrated
  path 'atac_integrated.rds', optional: true, emit: integrated
  script:
  """
  Rscript ${projectDir}/scripts/06b_preprocess_atac.R '${json}' ${objects.join(' ')}
  """
}

process MGATK_METRICS {
  tag "$id"
  cpus 2
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/mgatk", mode: 'copy'
  input:
  tuple val(id), val(mgatk), path(qc)
  output:
  path "${id}.mtDNA.tsv"
  script:
  """
  Rscript ${projectDir}/scripts/06c_mgatk_metrics.R '${id}' '${qc}' '${mgatk}' '${id}.mtDNA.tsv'
  """
}

process WNN_STAGE {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy'
  input:
  path rna
  path atac
  path json
  output:
  path 'wnn.rds'
  script:
  """
  Rscript ${projectDir}/scripts/08_wnn_stage.R wnn.rds '${json}' '${rna}' '${atac}'
  """
}

process ADD_MGATK {
  cpus 2
  memory '64 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy'
  input:
  path wnn
  path metrics
  output:
  path 'wnn_mgatk.rds'
  path 'all_cells_mtDNA.tsv'
  script:
  """
  Rscript ${projectDir}/scripts/09_add_mgatk.R '${wnn}' wnn_mgatk.rds all_cells_mtDNA.tsv ${metrics.join(' ')}
  """
}

process DEG_STAGE {
  cpus 4
  memory '64 GB'
  publishDir "${params.outdir}/${params.stage}/deg", mode: 'copy'
  input:
  path obj
  output:
  path 'cluster_markers.tsv'
  script:
  """
  Rscript ${projectDir}/scripts/10_cluster_markers.R '${obj}'
  """
}

workflow {
  rows = Channel.fromPath(params.samples, checkIfExists: true)
    .splitCsv(header: true)
    .filter { it.stage?.trim() == params.stage }
    .ifEmpty { error "No rows for stage ${params.stage} in ${params.samples}" }
    .map { row ->
      for (key in ['sample_id','arc_out','mgatk_dir'])
        if (!row[key]?.trim()) error "Missing ${key} in manifest row: ${row}"
      def id = row.sample_id.trim()
      if (!(id ==~ /5019_[A-Za-z0-9_.-]+/)) error "Expected m.5019 sample_id starting 5019_: ${id}"
      def arc = file(row.arc_out.trim()).toAbsolutePath().toString()
      def mt = file(row.mgatk_dir.trim()).toAbsolutePath().toString()
      if (!file("${arc}/filtered_feature_bc_matrix.h5").exists() || !file("${arc}/atac_fragments.tsv.gz").exists())
        error "Missing ARC H5 or fragments: ${arc}"
      if (!file("${arc}/atac_fragments.tsv.gz.tbi").exists()) error "Missing fragment index: ${arc}"
      if (!java.nio.file.Files.isDirectory(file(mt))) error "Missing mgatk directory: ${mt}"
      tuple(id, arc, mt)
    }

  PEAK_SAMPLE(rows)
  COMMON_PEAKS(PEAK_SAMPLE.out.map { id, bed -> bed }.collect())
  QC_METRICS(rows.combine(COMMON_PEAKS.out))
  QC_REVIEW(QC_METRICS.out.map { id, tab, plot -> tab }.collect())

  if (params.phase == 'cluster_review') {
    json_review = Channel.value(file(params.qc_json))
    PRE_DOUBLET_REVIEW(rows.combine(COMMON_PEAKS.out).combine(json_review))
  }

  if (params.phase == 'analysis') {
    json = Channel.value(file(params.qc_json))
    QC_FILTER_DOUBLET(rows.combine(COMMON_PEAKS.out).combine(json))
    filtered = QC_FILTER_DOUBLET.out
    MERGE_INTEGRATE_RNA(filtered.map { id, mt, rna, atac, report -> rna }.collect(), json)
    MERGE_INTEGRATE_ATAC(filtered.map { id, mt, rna, atac, report -> atac }.collect(), json)
    MGATK_METRICS(filtered.map { id, mt, rna, atac, report -> tuple(id, mt, rna) })
    WNN_STAGE(MERGE_INTEGRATE_RNA.out.stage, MERGE_INTEGRATE_ATAC.out.stage, json)
    ADD_MGATK(WNN_STAGE.out, MGATK_METRICS.out.collect())
    DEG_STAGE(ADD_MGATK.out[0])
  }
}
