nextflow.enable.dsl = 2

params.samples = null
params.stage = null
params.outdir = 'results'
params.macs3_path = null
params.phase = 'review'
// Nextflow loads these sections from -params-file; tasks receive scoped val maps.
params.peaks = null
params.note = null
params.thresholds = null
params.sample_overrides = null
params.review_QC_cluster = null
params.doublet = null
params.RNA_processing = null
params.ATAC_processing = null
params.WNN = null
params.copy_number = null

// Quote every shell argument, including paths containing spaces or apostrophes.
def shellQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
// Serialize only this task's val settings. The supplied pipeline JSON is never staged.
def taskSettingsFile(settings) {
  "cat > process_settings.json <<'PROCESS_SETTINGS_JSON'\n" +
    groovy.json.JsonOutput.toJson(settings) + "\nPROCESS_SETTINGS_JSON"
}
// Track shared files as path inputs, while R keeps their original absolute paths.
// Staging creates symlinks; BAMs/fragments are not copied into task directories.
def arcDataFiles(arc, boolean summary = false) {
  def names = ['filtered_feature_bc_matrix.h5', 'atac_fragments.tsv.gz', 'atac_fragments.tsv.gz.tbi']
  if (summary) names << 'summary.csv'
  names.collect { file("${arc}/${it}", checkIfExists: true) }
}
def mgatkDataFiles(directory, boolean depthOnly = false) {
  // Skip macOS "._" metadata files, which would otherwise match the patterns below.
  def names = file(directory).toFile().listFiles()?.findAll { it.isFile() && !it.name.startsWith('._') } ?: []
  def patterns = depthOnly ? [/[.]depthTable[.]txt$/] :
    [/[.]A[.]txt[.]gz$/, /[.]C[.]txt[.]gz$/, /[.]G[.]txt[.]gz$/, /[.]T[.]txt[.]gz$/,
     /_refAllele[.]txt(?:[.]gz)?$/, /[.]depthTable[.]txt$/]
  patterns.collect { pattern ->
    def matches = names.findAll { it.name =~ pattern }
    if (matches.size() != 1)
      error "Expected exactly one mgatk file matching ${pattern} in ${directory}"
    file(matches[0].toPath())
  }.sort { it.name }
}
def copyNumberDataFiles(arc, enriched, original, settings) {
  def paths = [file("${arc}/atac_possorted_bam.bam", checkIfExists: true)]
  paths.addAll(mgatkDataFiles(enriched, true))
  if (original) paths.addAll(mgatkDataFiles(original, true))
  def map = file("${arc}/per_barcode_metrics.csv")
  if (map.exists()) paths << map
  def copy = settings.copy_number
  if (copy.reference_bed) paths << file(copy.reference_bed, checkIfExists: true)
  else {
    paths << file("${arc}/atac_fragments.tsv.gz", checkIfExists: true)
    if (copy.blacklist) paths << file(copy.blacklist, checkIfExists: true)
  }
  paths.unique().sort { it.toString() }
}
// Compute once in the workflow and pass as val: calling this function only
// inside script interpolation does not track its returned digest in the cache.
def scriptsHash() {
  def sourceDigest = java.security.MessageDigest.getInstance('SHA-256')
  file("${projectDir}/scripts").toFile().listFiles()
    .findAll { it.name.endsWith('.R') }.sort { it.name }.each {
      sourceDigest.update(it.name.getBytes('UTF-8'))
      sourceDigest.update(it.bytes)
    }
  return sourceDigest.digest().encodeHex().toString()

}

process PEAK_SAMPLE {
  tag "$id"
  cpus 4
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/peaks/sample", mode: 'copy'
  input:
  val rSourceHash
  tuple val(id), val(arc), path(arc_inputs, stageAs: 'arc_inputs/*')
  output:
  tuple val(id), path("${id}.bed")
  script:
  """
  # R source checksum: ${rSourceHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/01_peaks.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(params.macs3_path)} '${id}.bed'
  """
}

process COMMON_PEAKS {
  cpus 2
  memory '12 GB'
  publishDir "${params.outdir}/${params.stage}/peaks", mode: 'copy'
  input:
  val rSourceHash
  path beds
  val peakMin
  val peakMax
  output:
  path 'common_peaks.bed'
  script:
  """
  # R source checksum: ${rSourceHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/02_consensus_peaks.R")} common_peaks.bed ${peakMin} ${peakMax} ${beds.collect { shellQuote(it) }.join(' ')}
  """
}

process QC_METRICS {
  tag "$id"
  cpus 4
  memory '48 GB'
  input:
  val rSourceHash
  tuple val(id), val(arc), path(peaks), path(arc_inputs, stageAs: 'arc_inputs/*')
  output:
  tuple val(id), path("${id}.metrics.tsv")
  script:
  """
  # R source checksum: ${rSourceHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/03_qc_metrics.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} '${id}.metrics.tsv'
  """
}

process QC_REVIEW {
  cpus 2
  memory '12 GB'
  publishDir "${params.outdir}/${params.stage}/qc/review", mode: 'copy'
  input:
  val rSourceHash
  path metric_files
  output:
  path 'qc_summary.tsv'
  path 'qc_mad_candidates.tsv'
  path 'all_cells_qc_metrics.tsv'
  path 'stage_qc_distributions.pdf'
  script:
  """
  # R source checksum: ${rSourceHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/04_qc_review.R")} qc_summary.tsv qc_mad_candidates.tsv stage_qc_distributions.pdf ${metric_files.collect { shellQuote(it) }.join(' ')}
  """
}

process PRE_DOUBLET_REVIEW {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/cluster_review", mode: 'copy'
  input:
  val rSourceHash
  tuple val(id), val(arc), path(peaks), val(settings), path(arc_inputs, stageAs: 'arc_inputs/*')
  output:
  path "${id}.cluster_cells.tsv"
  path "${id}.cluster_summary.tsv"
  path "${id}.cluster_qc.pdf"
  path "${id}.cluster_markers.tsv"
  script:
  """
  # R source checksum: ${rSourceHash}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/05_review_clusters.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} process_settings.json
  """
}

process QC_FILTER_DOUBLET {
  tag "$id"
  cpus 4
  memory '48 GB'
  publishDir "${params.outdir}/${params.stage}/qc/filtered", mode: 'copy', pattern: '*.{qc.tsv,pdf}'
  input:
  val rSourceHash
  tuple val(id), val(arc), path(peaks), val(settings), path(review), path(arc_inputs, stageAs: 'arc_inputs/*')
  output:
  tuple val(id), path("${id}.rna.qc.rds"), path("${id}.atac.qc.rds"), path("${id}.qc.tsv"), path("${id}_Doublet.pdf"), emit: objects
  tuple val(id), path("${id}.singlet_barcodes.tsv"), emit: barcodes
  script:
  """
  # R source checksum: ${rSourceHash}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/06_qc_doubletfinder.R")} ${shellQuote(id)} ${shellQuote(arc)} ${shellQuote(peaks)} process_settings.json '${id}.rna.qc.rds' '${id}.atac.qc.rds' '${id}.qc.tsv' ${shellQuote(review)}
  """
}

process MERGE_INTEGRATE_RNA {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: 'rna_*'
  input:
  val rSourceHash
  path objects
  val settings
  output:
  path 'rna_cell_annotations.tsv', emit: annotations
  path 'rna_unintegrated.rds', emit: unintegrated
  path 'rna_integrated.rds', optional: true, emit: integrated
  script:
  """
  # R source checksum: ${rSourceHash}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/07_preprocess_rna.R")} process_settings.json ${objects.collect { shellQuote(it) }.join(' ')}
  """
}

process MERGE_INTEGRATE_ATAC {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: 'atac_*.rds'
  input:
  val rSourceHash
  path objects
  val settings
  output:
  path 'atac_unintegrated.rds', emit: unintegrated
  path 'atac_integrated.rds', optional: true, emit: integrated
  script:
  """
  # R source checksum: ${rSourceHash}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/08_preprocess_atac.R")} process_settings.json ${objects.collect { shellQuote(it) }.join(' ')}
  """
}

process MGATK_METRICS {
  tag "$id"
  cpus 2
  memory '24 GB'
  publishDir "${params.outdir}/${params.stage}/mgatk", mode: 'copy'
  input:
  val rSourceHash
  tuple val(id), val(mgatk_enriched), path(qc), path(mgatk_inputs, stageAs: 'mgatk_inputs/*')
  output:
  tuple val(id), path("${id}.mtDNA.tsv")
  script:
  """
  # R source checksum: ${rSourceHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/09_mgatk_metrics.R")} ${shellQuote(id)} ${shellQuote(qc)} ${shellQuote(mgatk_enriched)} '${id}.mtDNA.tsv'
  """
}

// Step 10 uses step 06 singlets, independently of step 09 metric extraction.
process MTDNA_COPY_NUMBER {
  tag "$id"
  cpus 8
  memory '24 GB'
  time '48h'
  publishDir "${params.outdir}/${params.stage}/copy_number", mode: 'copy'
  input:
  val rSourceHash
  tuple val(id), val(arc), val(enriched), val(original), path(singlet_barcodes), val(settings), path(data_inputs, stageAs: 'external_inputs/input??/*')
  output:
  tuple val(id), path("${id}.mtDNA_copy_number.tsv"), emit: metrics
  tuple val(id), path("${id}.nuclear_reference.bed"), emit: reference
  tuple val(id), path("${id}.settings.txt"), emit: settings
  script:
  def ordinaryArg = original ? "--mgatk-dir " + shellQuote(original) : ''
  def mapFile = "${arc}/per_barcode_metrics.csv"
  def mapArg = file(mapFile).exists() ? "--barcode-map " + shellQuote(mapFile) : ''
  """
  # R source checksum: ${rSourceHash}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/10_estimate_mtdna_copy_number.R")} --sample ${shellQuote(id)} --bam ${shellQuote(arc + "/atac_possorted_bam.bam")} --fragments ${shellQuote(arc + "/atac_fragments.tsv.gz")} --metrics-tsv ${shellQuote(singlet_barcodes)} --mgatk-enriched-dir ${shellQuote(enriched)} ${ordinaryArg} ${mapArg} --json process_settings.json --out-prefix ${shellQuote(id)}
  """
}

process WNN_STAGE {
  cpus 8
  memory '96 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy', pattern: '*.rds'
  publishDir "${params.outdir}/${params.stage}/plots", mode: 'copy', pattern: '*.pdf'
  input:
  val rSourceHash
  path rna
  path atac
  val settings
  output:
  path 'wnn.rds', emit: object
  path '*_umap.pdf', emit: plots
  script:
  """
  # R source checksum: ${rSourceHash}
  ${taskSettingsFile(settings)}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/12_wnn_stage.R")} wnn.rds process_settings.json ${shellQuote(rna)} ${shellQuote(atac)}
  """
}

process ADD_MGATK {
  cpus 2
  memory '64 GB'
  publishDir "${params.outdir}/${params.stage}/objects", mode: 'copy'
  input:
  val rSourceHash
  path wnn
  path metrics
  path copy_numbers
  output:
  path 'wnn_mgatk.rds'
  path 'all_cells_mtDNA.tsv'
  script:
  """
  # R source checksum: ${rSourceHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/13_add_mgatk.R")} ${shellQuote(wnn)} wnn_mgatk.rds all_cells_mtDNA.tsv ${metrics.collect { shellQuote(it) }.join(' ')} ${copy_numbers ? "--copy-number-files " + copy_numbers.collect { shellQuote(it) }.join(' ') : ''}
  """
}

process DEG_STAGE {
  cpus 4
  memory '64 GB'
  publishDir "${params.outdir}/${params.stage}/deg", mode: 'copy'
  input:
  val rSourceHash
  path obj
  output:
  path 'cluster_markers.tsv'
  script:
  """
  # R source checksum: ${rSourceHash}
  Rscript ${shellQuote(projectDir.toString() + "/scripts/14_cluster_markers.R")} ${shellQuote(obj)}
  """
}

workflow {
  def rSourceHash = scriptsHash()
  if (!params.samples || !params.stage) error 'Set samples and stage in the file supplied with -params-file'
  if (!(params.phase in ['review', 'cluster_review', 'analysis'])) error '--phase must be review, cluster_review or analysis'
  if (!(params.stage.toString() ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'Stage must be a safe directory name'
  // Nextflow parses -params-file. Keep run paths out of the analysis task maps.
  // Initial review reads peaks while still permitting unfinished QC bounds.
  def pipelineSettings = PipelineSettings.analysis(new LinkedHashMap(params))
  // Peak widths come from the pipeline JSON; omitted settings use defaults.
  def selectedPeaks = PipelineSettings.peaks(pipelineSettings.peaks)
  def peakMin = selectedPeaks.peak_min
  def peakMax = selectedPeaks.peak_max
  if (!params.macs3_path) error 'Set macs3_path to the absolute MACS3 executable path in -params-file'
  if (!java.nio.file.Paths.get(params.macs3_path.toString()).isAbsolute()) error 'macs3_path must be absolute'
  if (!java.nio.file.Files.isRegularFile(file(params.macs3_path)) ||
      !java.nio.file.Files.isExecutable(file(params.macs3_path)))
    error "MACS3 executable not found or not executable: ${params.macs3_path}"
  if (params.enable_conda) {
    if (!params.conda_env) error 'The conda profile requires conda_env in -params-file (shared environment path or YAML)'
    if (!file(params.conda_env).exists()) error "Conda environment or YAML not found: ${params.conda_env}"
  }

  // Validate the complete supplied config once before projecting task inputs.
  // Initial review intentionally permits unfinished QC bounds.
  def scopedSettings = params.phase == 'review' ? [:] : PipelineSettings.normalize(pipelineSettings)
  def rnaSettings = params.phase == 'analysis' ? PipelineSettings.rna(scopedSettings) : [:]
  def atacSettings = params.phase == 'analysis' ? PipelineSettings.atac(scopedSettings) : [:]
  def wnnSettings = params.phase == 'analysis' ? PipelineSettings.wnn(scopedSettings) : [:]
  def copyNumberSettings = params.phase == 'analysis' ? PipelineSettings.copyNumber(scopedSettings) : [:]

  // Stop on misspelled sample IDs, which would otherwise be silently ignored.
  if (params.phase != 'review') {
    def stageIds = file(params.samples).splitCsv(header: true)
      .findAll { it.stage?.trim() == params.stage }.collect { it.sample_id?.trim() }
    def jsonIds = scopedSettings.sample_overrides.keySet() +
      scopedSettings.doublet.sample_rates.keySet() + scopedSettings.doublet.exclude_clusters.keySet()
    def unknown = jsonIds.findAll { !(it in stageIds) }.unique()
    if (unknown) error "Sample IDs in pipeline JSON not found for stage ${params.stage}: ${unknown.join(', ')}"
  }

  seenIds = new HashSet()
  // Preserve both mitochondrial sources; step 09 uses enriched, step 10 uses both.
  sample_records = Channel.fromPath(params.samples, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
      if (row.keySet() != ['stage', 'sample_id', 'arc_out', 'mgatk_enriched_dir', 'mgatk_dir'].toSet())
        error 'Sample CSV must contain exactly stage,sample_id,arc_out,mgatk_enriched_dir,mgatk_dir columns'
      row
    }
    .filter { it.stage?.trim() == params.stage }
    .ifEmpty { error "No rows for stage ${params.stage} in ${params.samples}" }
    .map { row ->
      ['sample_id', 'arc_out', 'mgatk_enriched_dir'].each { key ->
        if (!row[key]?.trim()) error "Missing ${key} in manifest row: ${row}"
      }
      def id = row.sample_id.trim()
      if (!(id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
        error "Invalid sample_id ${id}: start with a letter or number and use only letters, numbers, _, . or -"
      if (id == 'ALL_STAGE') error 'ALL_STAGE is reserved for the combined QC summary'
      if (!seenIds.add(id)) error "Duplicate sample_id for stage ${params.stage}: ${id}"
      def arc = file(row.arc_out.trim()).toAbsolutePath().toString()
      def enriched = file(row.mgatk_enriched_dir.trim()).toAbsolutePath().toString()
      // Blank ordinary mgatk is allowed: its estimate will be NA until supplied.
      def original = row.mgatk_dir?.trim() ? file(row.mgatk_dir.trim()).toAbsolutePath().toString() : ''
      if (!file("${arc}/filtered_feature_bc_matrix.h5").exists() || !file("${arc}/atac_fragments.tsv.gz").exists())
        error "Missing ARC H5 or fragments: ${arc}"
      if (!file("${arc}/atac_fragments.tsv.gz.tbi").exists()) error "Missing fragment index: ${arc}"
      if (!java.nio.file.Files.isDirectory(file(enriched)))
        error "Missing mgatk_enriched_dir directory: ${enriched}"
      if (params.phase == 'analysis') {
        if (!file("${arc}/summary.csv").exists()) error "Missing ARC summary.csv: ${arc}"
        mgatkDataFiles(enriched)
        if (!file("${arc}/atac_possorted_bam.bam").exists())
          error "Step 10 requires ARC atac_possorted_bam.bam: ${arc}"
        if (original && !java.nio.file.Files.isDirectory(file(original)))
          error "Missing mgatk_dir directory: ${original}"
        copyNumberDataFiles(arc, enriched, original, copyNumberSettings)
      }
      tuple(id, arc, enriched, original)
    }
    .toSortedList { a, b -> a[0] <=> b[0] }
    .flatMap { it }

  rows = sample_records.map { id, arc, enriched, original -> tuple(id, arc, enriched) }

  peak_inputs = rows.map { id, arc, enriched -> tuple(id, arc, arcDataFiles(arc)) }
  PEAK_SAMPLE(rSourceHash, peak_inputs)
  COMMON_PEAKS(rSourceHash, PEAK_SAMPLE.out.map { id, bed -> bed }.collect().map { it.sort { path -> path.name } }, peakMin, peakMax)
  metric_inputs = rows.combine(COMMON_PEAKS.out).map { id, arc, enriched, peaks ->
    tuple(id, arc, peaks, arcDataFiles(arc))
  }
  QC_METRICS(rSourceHash, metric_inputs)
  QC_REVIEW(rSourceHash, QC_METRICS.out.map { id, tab -> tab }.collect().map { it.sort { path -> path.name } })

  if (params.phase == 'cluster_review') {
    review_inputs = rows.combine(COMMON_PEAKS.out).map { id, arc, enriched, peaks ->
      tuple(id, arc, peaks, PipelineSettings.sample(scopedSettings, id), arcDataFiles(arc))
    }
    PRE_DOUBLET_REVIEW(rSourceHash, review_inputs)
  }

  if (params.phase == 'analysis') {
    def reviewDir = "${params.outdir}/${params.stage}/qc/cluster_review"
    reviewed = rows.combine(COMMON_PEAKS.out).map { id, arc, enriched, peaks ->
      def reviewFile = file("${reviewDir}/${id}.cluster_cells.tsv", checkIfExists: true)
      tuple(id, arc, peaks, PipelineSettings.sample(scopedSettings, id, true), reviewFile, arcDataFiles(arc, true))
    }
    QC_FILTER_DOUBLET(rSourceHash, reviewed)
    filtered = QC_FILTER_DOUBLET.out.objects
    MERGE_INTEGRATE_RNA(rSourceHash, filtered.map { id, rna, atac, report, plot -> rna }.collect().map { it.sort { path -> path.name } }, rnaSettings)
    MERGE_INTEGRATE_ATAC(rSourceHash, filtered.map { id, rna, atac, report, plot -> atac }.collect().map { it.sort { path -> path.name } }, atacSettings)
    mgatk_inputs = sample_records.join(filtered).map { id, arc, enriched, original, rna, atac, report, plot ->
      tuple(id, enriched, rna, mgatkDataFiles(enriched))
    }
    MGATK_METRICS(rSourceHash, mgatk_inputs)
    // Start copy-number tasks from each sample's singlets without waiting for
    // RNA/ATAC integration or the heteroplasmy tables produced by step 09.
    copy_inputs = sample_records.join(QC_FILTER_DOUBLET.out.barcodes).map { id, arc, enriched, original, barcodes ->
      tuple(id, arc, enriched, original, barcodes, copyNumberSettings,
        copyNumberDataFiles(arc, enriched, original, copyNumberSettings))
    }
    MTDNA_COPY_NUMBER(rSourceHash, copy_inputs)
    copy_numbers = MTDNA_COPY_NUMBER.out.metrics.map { id, tab -> tab }.collect().map { it.sort { path -> path.name } }
    // Use the saved analysis object directly; omitted switches default to true.
    rna_for_wnn = scopedSettings.RNA_processing.integrate_rna == false ?
      MERGE_INTEGRATE_RNA.out.unintegrated : MERGE_INTEGRATE_RNA.out.integrated
    atac_for_wnn = scopedSettings.ATAC_processing.integrate_atac == false ?
      MERGE_INTEGRATE_ATAC.out.unintegrated : MERGE_INTEGRATE_ATAC.out.integrated
    WNN_STAGE(rSourceHash, rna_for_wnn, atac_for_wnn, wnnSettings)
    ADD_MGATK(rSourceHash, WNN_STAGE.out.object, MGATK_METRICS.out.map { id, tab -> tab }.collect().map { it.sort { path -> path.name } }, copy_numbers)
    // RNA cluster markers need WNN identities, independently of mtDNA metadata.
    DEG_STAGE(rSourceHash, WNN_STAGE.out.object)
  }
}
