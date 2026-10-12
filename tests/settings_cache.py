#!/usr/bin/env python3
"""Exercise the real Nextflow graph/cache with mock R outputs, without sequencing data.

Run: JAVA_CMD=/path/to/java17 python3 tests/settings_cache.py --nextflow /path/to/nextflow
All inputs, outputs and caches stay in a temporary directory. This checks task
selection and resume behavior, not Seurat/Signac analysis or copy-number accuracy.
"""
import argparse
import csv
import json
import os
import shutil
from pathlib import Path
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[1]
MOCK = r'''#!/usr/bin/env python3
import json,os,sys,time
from pathlib import Path
a=sys.argv[2:]; step=Path(sys.argv[1]).name[:2]
# In the baseline run, keep integration and mgatk tasks open until the sibling
# copy-number tasks start/finish. An accidental dependency then fails promptly.
parallel=os.environ.get('STAGE_TEST_PARALLEL_DIR'); label=None
if parallel and step in ['07','08','09','10']:
 root=Path(parallel)
 label=step if step in ['07','08'] else step+'-'+(a[0] if step=='09' else a[a.index('--sample')+1])
 (root/(label+'.started')).touch()
 required=[f'10-{a[0]}.started'] if step=='09' else (
  [f'{s}-{id}.finished' for s in ['09','10'] for id in ['WT_1','WT_2']] if step in ['07','08'] else [])
 deadline=time.monotonic()+30
 while not all((root/name).exists() for name in required):
  if time.monotonic()>deadline: raise RuntimeError('Sibling mitochondrial tasks did not run independently: '+label)
  time.sleep(0.1)
settings=json.loads(Path('process_settings.json').read_text()) if Path('process_settings.json').exists() else {}
body=json.dumps(settings,sort_keys=True)+'\n'
def write(name): Path(name).write_text(body)
if step=='01': write(a[3])
elif step=='02': write(a[0])
elif step=='03': write(a[3])
elif step=='04':
 for name in a[:3]+['all_cells_qc_metrics.tsv']: write(name)
elif step=='05':
 for suffix in ['cluster_cells.tsv','cluster_summary.tsv','cluster_qc.pdf','cluster_markers.tsv']: write(a[0]+'.'+suffix)
elif step=='06':
 for name in a[4:7]+[a[0]+'_Doublet.pdf']: write(name)
 write(a[0]+'.singlet_barcodes.tsv')
elif step in ['07','08']:
 prefix='rna' if step=='07' else 'atac'; section='RNA_processing' if step=='07' else 'ATAC_processing'
 write(prefix+'_unintegrated.rds')
 if settings[section]['integrate_'+prefix]: write(prefix+'_integrated.rds')
 if step=='07': write('rna_cell_annotations.tsv')
elif step=='09': write(a[3])
elif step=='10':
 prefix=a[a.index('--out-prefix')+1]
 for suffix in ['mtDNA_copy_number.tsv','nuclear_reference.bed','settings.txt']: write(prefix+'.'+suffix)
elif step=='12':
 write(a[0]); write('wnn_clusters_umap.pdf'); write('Cell_lineage_umap.pdf')
elif step=='13': write(a[1]); write(a[2])
elif step=='14': write('cluster_markers.tsv')
else: raise RuntimeError('Unexpected script '+sys.argv[1])
if label: (root/(label+'.finished')).touch()
'''

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--nextflow', default='nextflow')
    ap.add_argument('--keep', action='store_true', help='Retain synthetic run directory for inspection')
    args = ap.parse_args()
    root = Path(tempfile.mkdtemp(prefix='stage-settings-cache-'))
    print(f'Synthetic runs: {root}', flush=True)
    try:
        mockbin = root / 'bin'; mockbin.mkdir()
        executable = mockbin / 'Rscript'; executable.write_text(MOCK); executable.chmod(0o755)
        macs = mockbin / 'macs3'; macs.write_text('#!/bin/sh\nexit 0\n'); macs.chmod(0o755)
        manifest = root / 'samples.csv'
        with manifest.open('w') as handle:
            writer = csv.writer(handle)
            writer.writerow(['stage', 'sample_id', 'arc_out', 'mgatk_enriched_dir', 'mgatk_dir'])
            for sample in ['WT_1', 'WT_2']:
                arc = root / sample / 'arc'; arc.mkdir(parents=True)
                for filename in ['filtered_feature_bc_matrix.h5', 'atac_fragments.tsv.gz',
                                 'atac_fragments.tsv.gz.tbi', 'atac_possorted_bam.bam', 'summary.csv']:
                    (arc / filename).write_text('synthetic\n')
                enriched = root / sample / 'enriched'; enriched.mkdir()
                for name in ['sample.A.txt.gz', 'sample.C.txt.gz', 'sample.G.txt.gz',
                             'sample.T.txt.gz', 'chrM_refAllele.txt', 'sample.depthTable.txt']:
                    (enriched / name).write_text('synthetic\n')
                original = root / sample / 'ordinary'; original.mkdir()
                (original / 'sample.depthTable.txt').write_text('synthetic\n')
                writer.writerow(['E8.5', sample, arc, enriched, original])
        config = root / 'test.config'
        config.write_text("process { withName: '.*' { cpus = 1; memory = '128 MB' } }\n"
                          "executor.queueSize = 8\nexecutor.cpus = 8\n")
        settings_file = root / 'settings.json'
        settings = json.loads((REPO / 'examples/pipeline_settings.example.json').read_text())
        settings.update(samples=str(manifest), stage='E8.5', outdir=str(root / 'results'),
                        macs3_path=str(macs), conda_env=None)
        env = dict(os.environ, PATH=str(mockbin) + os.pathsep + os.environ['PATH'], NXF_ANSI_LOG='false')
        def run(label, phase='analysis', parallel=False):
            settings_file.write_text(json.dumps(settings, indent=2))
            command = [args.nextflow, '-log', str(root / (label + '.nextflow.log')), 'run', str(REPO / 'main.nf'),
                       '-c', str(config), '-work-dir', str(root / 'work'),
                       '-params-file', str(settings_file), '--phase', phase, '-resume']
            run_env = dict(env)
            if parallel:
                parallel_dir = root / 'parallel'; parallel_dir.mkdir()
                run_env['STAGE_TEST_PARALLEL_DIR'] = str(parallel_dir)
            result = subprocess.run(command, cwd=root, env=run_env, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=180)
            (root / (label + '.output.log')).write_text(result.stdout)
            if result.returncode:
                raise AssertionError(f'{label} failed:\n{result.stdout}')
            trace = root / 'results' / 'E8.5' / 'reports' / f'execution_trace_{phase}.txt'
            rows = list(csv.DictReader(trace.open(), delimiter='\t'))
            (root / (label + '.trace.tsv')).write_text(trace.read_text())
            completed = {row['name'] for row in rows if row['status'] == 'COMPLETED'}
            print(label + ': ran ' + ', '.join(sorted(completed)), flush=True)
            return completed
        def expect(actual, processes):
            # Untagged stages have their process name; per-sample names end in (ID).
            normalized = {name.split(' (')[0] for name in actual}
            assert normalized == set(processes), (actual, processes)
        run('initial-review', 'review')
        changed = run('cluster-review', 'cluster_review')
        assert changed == {'PRE_DOUBLET_REVIEW (WT_1)', 'PRE_DOUBLET_REVIEW (WT_2)'}, changed
        run('baseline', parallel=True)
        markers = root / 'parallel'
        for sample in ['WT_1', 'WT_2']:
            copy_started = (markers / f'10-{sample}.started').stat().st_mtime_ns
            mgatk_finished = (markers / f'09-{sample}.finished').stat().st_mtime_ns
            assert copy_started < mgatk_finished
            for stage in ['07', '08']:
                assert mgatk_finished < (markers / f'{stage}.finished').stat().st_mtime_ns
        print('Mgatk and copy-number tasks overlap while RNA/ATAC integration waits.', flush=True)
        expect(run('unchanged'), [])
        settings['WNN']['resolution'] = 0.8
        expect(run('wnn-resolution'), ['WNN_STAGE', 'ADD_MGATK', 'DEG_STAGE'])
        expect(run('cluster-review-wnn-edit', 'cluster_review'), [])
        # Already-computed PC50 covers this range; RNA preprocessing stays cached.
        settings['WNN']['RNA_dims'] = '1:30'
        expect(run('wnn-covered-dims'), ['WNN_STAGE', 'ADD_MGATK', 'DEG_STAGE'])
        settings['WNN']['RNA_dims'] = '1:55'
        expect(run('wnn-more-pcs'), ['MERGE_INTEGRATE_RNA', 'WNN_STAGE', 'ADD_MGATK', 'DEG_STAGE'])
        # JSON formatting/key order and note do not affect any task inputs.
        settings = dict(reversed(list(settings.items())))
        settings['note'] = 'This note does not control analysis.'
        expect(run('format-note'), [])
        settings['doublet']['sample_rates'] = {'WT_1': 0.04}
        changed = run('sample-doublet')
        expect(changed, ['QC_FILTER_DOUBLET', 'MERGE_INTEGRATE_RNA', 'MERGE_INTEGRATE_ATAC',
                         'MGATK_METRICS', 'MTDNA_COPY_NUMBER', 'WNN_STAGE', 'ADD_MGATK', 'DEG_STAGE'])
        assert 'QC_FILTER_DOUBLET (WT_1)' in changed and 'QC_FILTER_DOUBLET (WT_2)' not in changed
        expect(run('cluster-review-doublet-edit', 'cluster_review'), [])
        # Settings for a sample ID absent from the manifest are rejected, not ignored.
        settings['doublet']['sample_rates']['unselected'] = 0.07
        try:
            run('unselected-sample')
            raise RuntimeError('Unknown sample ID was accepted')
        except AssertionError as error:
            assert 'Sample IDs in pipeline JSON not found' in str(error), error
        del settings['doublet']['sample_rates']['unselected']
        settings['copy_number']['ploidy'] = 4
        expect(run('copy-number'), ['MTDNA_COPY_NUMBER', 'ADD_MGATK'])
        # Integration switches select actual outputs/reductions at workflow level.
        settings['ATAC_processing']['integrate_atac'] = False
        expect(run('disable-atac-integration'), ['MERGE_INTEGRATE_ATAC', 'WNN_STAGE', 'ADD_MGATK', 'DEG_STAGE'])
        settings['ATAC_processing']['ATAC_integrate_dims'] = '1:60'
        expect(run('unused-atac-dims'), [])
        settings['sample_overrides'] = {'WT_1': {'nCount_RNA': {'min': 1200, 'max': 40000}}}
        changed = run('sample-qc-review', 'cluster_review')
        assert changed == {'PRE_DOUBLET_REVIEW (WT_1)'}, changed
        changed = run('sample-qc-analysis')
        assert 'QC_FILTER_DOUBLET (WT_1)' in changed and 'QC_FILTER_DOUBLET (WT_2)' not in changed
        # Edits at the same shared path must invalidate the relevant file inputs.
        (root / 'WT_1' / 'ordinary' / 'sample.depthTable.txt').write_text('changed ordinary depth\n')
        changed = run('ordinary-depth-edit')
        expect(changed, ['MTDNA_COPY_NUMBER', 'ADD_MGATK'])
        assert 'MTDNA_COPY_NUMBER (WT_1)' in changed and 'MTDNA_COPY_NUMBER (WT_2)' not in changed
        (root / 'WT_1' / 'enriched' / 'sample.A.txt.gz').write_text('changed allele counts\n')
        changed = run('enriched-count-edit')
        expect(changed, ['MGATK_METRICS', 'ADD_MGATK'])
        assert 'MGATK_METRICS (WT_1)' in changed and 'MGATK_METRICS (WT_2)' not in changed
        (root / 'WT_1' / 'arc' / 'summary.csv').write_text('changed detected cells\n')
        changed = run('arc-summary-edit')
        expect(changed, ['QC_FILTER_DOUBLET', 'MERGE_INTEGRATE_RNA', 'MERGE_INTEGRATE_ATAC',
                         'MGATK_METRICS', 'MTDNA_COPY_NUMBER', 'WNN_STAGE', 'ADD_MGATK', 'DEG_STAGE'])
        assert 'QC_FILTER_DOUBLET (WT_1)' in changed and 'QC_FILTER_DOUBLET (WT_2)' not in changed
        (root / 'WT_1' / 'arc' / 'atac_possorted_bam.bam').write_text('changed BAM input\n')
        changed = run('arc-bam-edit')
        expect(changed, ['MTDNA_COPY_NUMBER', 'ADD_MGATK'])
        assert 'MTDNA_COPY_NUMBER (WT_1)' in changed and 'MTDNA_COPY_NUMBER (WT_2)' not in changed
        reference = root / 'shared-reference.bed'; reference.write_text('chr1\t0\t100\n')
        settings['copy_number']['reference_bed'] = str(reference)
        expect(run('shared-reference'), ['MTDNA_COPY_NUMBER', 'ADD_MGATK'])
        reference.write_text('chr1\t0\t1200\n')
        expect(run('shared-reference-edit'), ['MTDNA_COPY_NUMBER', 'ADD_MGATK'])
        settings['copy_number']['quantile'] = 0.9
        expect(run('unused-reference-quantile'), [])
        # Optional ordinary mgatk changes only the sample's copy-number task.
        with manifest.open() as handle:
            records = list(csv.DictReader(handle))
        for record in records:
            if record['sample_id'] == 'WT_2': record['mgatk_dir'] = ''
        with manifest.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader(); writer.writerows(records)
        changed = run('blank-ordinary-source')
        expect(changed, ['MTDNA_COPY_NUMBER', 'ADD_MGATK'])
        assert 'MTDNA_COPY_NUMBER (WT_2)' in changed and 'MTDNA_COPY_NUMBER (WT_1)' not in changed
        # Changed ARC reads require new peaks/common intervals and cluster review.
        (root / 'WT_1' / 'arc' / 'filtered_feature_bc_matrix.h5').write_text('changed H5 input\n')
        changed = run('arc-h5-review', 'cluster_review')
        expect(changed, ['PEAK_SAMPLE', 'COMMON_PEAKS', 'QC_METRICS', 'QC_REVIEW', 'PRE_DOUBLET_REVIEW'])
        assert 'PEAK_SAMPLE (WT_1)' in changed and 'PEAK_SAMPLE (WT_2)' not in changed
        # Check that the real R reader accepts settings generated by this workflow.
        real_r = shutil.which('Rscript')
        if real_r:
            subprocess.run([real_r, str(REPO / 'tests/process_settings.R'), str(root / 'work')],
                           cwd=REPO, check=True, timeout=30)
        print('All Nextflow settings/cache checks passed.', flush=True)
    finally:
        if not args.keep:
            shutil.rmtree(root)

if __name__ == '__main__':
    main()
