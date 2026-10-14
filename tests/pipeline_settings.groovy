// Run with Groovy (the Nextflow distribution includes Groovy).
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

def PipelineSettings = new GroovyClassLoader(getClass().classLoader).parseClass(new File('lib/PipelineSettings.groovy'))

def original = new JsonSlurper().parse(new File('examples/pipeline_settings.example.json'))
def fresh = { new JsonSlurper().parseText(JsonOutput.toJson(original)) }
def expectError = { Closure action ->
    boolean failed = false
    try { action() } catch (IllegalArgumentException ignored) { failed = true }
    assert failed
}
assert !PipelineSettings.analysis(original).keySet().any { it in PipelineSettings.RUN_PARAMETERS }
expectError { PipelineSettings.analysis(original + [json: '/old/settings.json']) }
expectError { PipelineSettings.analysis(original + [cluster_review_dir: '/old/review']) }
def changedPaths = fresh(); changedPaths.outdir = '/different/results'; changedPaths.phase = 'analysis'
assert PipelineSettings.normalize(changedPaths) == PipelineSettings.normalize(original)
def base = PipelineSettings.normalize(fresh())
assert PipelineSettings.rna(base).component_count == 50
assert PipelineSettings.atac(base).component_count == 40
assert PipelineSettings.wnn(base).RNA_processing.keySet() == ['integrate_rna'].toSet()
assert PipelineSettings.wnn(base).ATAC_processing.keySet() == ['integrate_atac'].toSet()
def modified = fresh()
modified.WNN.resolution = 0.5
modified.WNN.RNA_dims = '1:30'
modified.note = 'Only for humans'
def changed = PipelineSettings.normalize(modified)
assert PipelineSettings.sample(base, 'WT_1', true) == PipelineSettings.sample(changed, 'WT_1', true)
assert PipelineSettings.rna(base) == PipelineSettings.rna(changed)
assert PipelineSettings.atac(base) == PipelineSettings.atac(changed)
assert PipelineSettings.copyNumber(base) == PipelineSettings.copyNumber(changed)
modified.WNN.RNA_dims = '1:60'
assert PipelineSettings.rna(PipelineSettings.normalize(modified)).component_count == 60
// Resolve only effective current-sample bounds, rates and cluster exclusions.
modified = fresh()
modified.sample_overrides = [WT_1: [nCount_RNA: [min: 1500, max: null]]]
modified.doublet.expected_rate = 0.08
modified.doublet.sample_rates = [WT_1: 0.04, WT_2: 'auto']
modified.doublet.exclude_clusters = [WT_1: ['2']]
changed = PipelineSettings.normalize(modified)
assert PipelineSettings.sample(changed, 'WT_1', true).thresholds.nCount_RNA == [min: 1500, max: null]
assert PipelineSettings.sample(changed, 'WT_2', true).thresholds == base.thresholds
assert PipelineSettings.sample(changed, 'WT_1', true).doublet.expected_rate == 0.04
assert PipelineSettings.sample(changed, 'WT_2', true).doublet.expected_rate == 'auto'
assert PipelineSettings.sample(changed, 'WT_1', true).doublet.exclude_clusters == [WT_1: ['2']]
modified.doublet.expected_rate = 0.06
assert PipelineSettings.sample(changed, 'WT_1', true) == PipelineSettings.sample(PipelineSettings.normalize(modified), 'WT_1', true)
assert PipelineSettings.sample(changed, 'WT_1') == PipelineSettings.sample(PipelineSettings.normalize(modified), 'WT_1')
// Settings inactive in the chosen mode never enter task maps.
modified = fresh(); modified.ATAC_processing.integrate_atac = false
changed = PipelineSettings.normalize(modified)
assert !PipelineSettings.atac(changed).ATAC_processing.containsKey('ATAC_integrate_dims')
modified.copy_number.reference_bed = '/shared/reference.bed'
changed = PipelineSettings.normalize(modified)
assert !PipelineSettings.copyNumber(changed).copy_number.containsKey('bin_size')
assert !PipelineSettings.copyNumber(changed).copy_number.containsKey('quantile')
assert !PipelineSettings.copyNumber(changed).copy_number.containsKey('blacklist')
// Integer-valued JSON widths and null PCs agree with direct R config defaults.
assert PipelineSettings.peaks([peak_min: 20.0, peak_max: 10000.0]) == [peak_min: 20, peak_max: 10000]
assert PipelineSettings.peaks([peak_min: 50, peak_max: 5000]) == [peak_min: 50, peak_max: 5000]
assert PipelineSettings.peaks(null) == [peak_min: 20, peak_max: 10000]
assert PipelineSettings.peaks([peak_min: 50]) == [peak_min: 50, peak_max: 10000]
expectError { PipelineSettings.peaks([peak_min: 20.5]) }
expectError { PipelineSettings.peaks([peak_min: '20']) }
expectError { PipelineSettings.peaks([peak_max: 2147483648L]) }
modified = fresh(); modified.doublet.PC_dims = null
assert PipelineSettings.sample(PipelineSettings.normalize(modified), 'WT_1', true).doublet.PC_dims == '1:20'
// remove_mf_doublets defaults to true and reaches each sample's doublet task.
modified = fresh(); modified.doublet.remove("remove_mf_doublets")
assert PipelineSettings.sample(PipelineSettings.normalize(modified), 'WT_1', true).doublet.remove_mf_doublets == true
modified = fresh(); modified.doublet.remove_mf_doublets = false
assert PipelineSettings.sample(PipelineSettings.normalize(modified), 'WT_1', true).doublet.remove_mf_doublets == false
// MGATK variant thresholds default to 0.01 / 0.65 / 300 and reach step 09's settings.
modified = fresh(); modified.remove('MGATK')
assert PipelineSettings.mgatk(PipelineSettings.normalize(modified)).MGATK ==
    [n_cells_conf_detected: 300, strand_correlation: 0.65, vmr: 0.01]
modified = fresh(); modified.MGATK.vmr = 0.05
assert PipelineSettings.mgatk(PipelineSettings.normalize(modified)).MGATK.vmr == 0.05
modified = fresh(); modified.sample_overrides = ['': [:]]
expectError { PipelineSettings.normalize(modified) }
// Fail before scheduling tasks, including invalid overrides for other samples.
for (mutate in [
    { it.thresholds.nCount_RNA = [min: 10, max: 2] },
    { it.thresholds.'percent.mt' = [min: null, max: null] },
    { it.sample_overrides = [other: [TSS: [min: 2, max: null]]] },
    { it.doublet.sample_rates = [other: 2] },
    { it.doublet.exclude_clusters = [other: [2]] },
    { it.doublet.PC_dims = '2:20' },
    { it.doublet.remove_mf_doublets = 'yes' },
    { it.MGATK.vmr = 0 },
    { it.MGATK.strand_correlation = 1.5 },
    { it.MGATK.n_cells_conf_detected = 2.5 },
    { it.MGATK.unknown = 1 },
    { it.RNA_processing.algorithm = '3' },
    { it.RNA_processing.resolution = 0 },
    { it.ATAC_processing.integrate_atac = 'false' },
    { it.ATAC_processing.ATAC_integrate_dims = '2:40' },
    { it.WNN.ATAC_dims = '2:50' },
    { it.WNN.RNA_dims = '1:501' },
    { it.WNN.plot_point_size = 0 },
    { it.WNN.lineage_markers = ['Pou3f1', 'Pou3f1'] },
    { it.WNN.lineage_assay = '' },
    { it.ATAC_processing.min_cutoff = 'q100' },
    { it.copy_number.mapq = 1.5 },
    { it.copy_number.reference_bed = 'relative.bed' },
    { it.copy_number.keep_duplicates = 'false' },
    { it.copy_number.quantile = 1 },
    { it.copy_number.ploidy = 0 },
    { it.Wnn = [:] },
    { it.RNA_processing.typo = 1 }
]) {
    modified = fresh(); mutate(modified)
    expectError { PipelineSettings.normalize(modified) }
}
// Omitting sections/keys uses the same defaults as explicit values.
modified = fresh()
['review_QC_cluster', 'doublet', 'RNA_processing', 'ATAC_processing', 'WNN', 'copy_number'].each { modified.remove(it) }
changed = PipelineSettings.normalize(modified)
assert PipelineSettings.sample(base, 'WT_1', true) == PipelineSettings.sample(changed, 'WT_1', true)
assert PipelineSettings.rna(base) == PipelineSettings.rna(changed)
assert PipelineSettings.atac(base) == PipelineSettings.atac(changed)
assert PipelineSettings.wnn(base) == PipelineSettings.wnn(changed)
assert PipelineSettings.copyNumber(base) == PipelineSettings.copyNumber(changed)
println 'All workflow settings projection checks passed.'
