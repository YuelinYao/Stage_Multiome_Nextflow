#!/usr/bin/env python3
"""Check R-source cache invalidation with real Nextflow in an isolated fixture.

Run: JAVA_CMD=/path/to/java17 python3 tests/script_hash.py --nextflow /path/to/nextflow
"""
import argparse
import csv
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow', default='nextflow')
    args = parser.parse_args()
    source = (Path(__file__).resolve().parents[1] / 'main.nf').read_text()
    helper = source[source.index('def scriptsHash() {'):source.index('\nprocess PEAK_SAMPLE')]
    assert source.count('  val rSourceHash\n') == 13
    assert source.count('# R source checksum: ${rSourceHash}') == 13
    with tempfile.TemporaryDirectory(prefix='stage-script-hash-') as directory:
        root = Path(directory)
        (root / 'scripts').mkdir()
        script = root / 'scripts' / 'example.R'
        script.write_text('message("original")\n')
        (root / 'main.nf').write_text('nextflow.enable.dsl=2\n' + helper + '''
process CHECK_HASH {
  input:
  val rSourceHash
  output:
  path 'result.txt'
  script:
  """
  # R source checksum: ${rSourceHash}
  echo checked > result.txt
  """
}
workflow {
  def rSourceHash = scriptsHash()
  CHECK_HASH(rSourceHash)
}
''')
        hashes = []
        for label, expected in [('initial', 'COMPLETED'), ('unchanged', 'CACHED'),
                                ('edited', 'COMPLETED')]:
            if label == 'edited':
                script.write_text('message("edited")\n')
            trace = root / (label + '.trace.tsv')
            result = subprocess.run([args.nextflow, 'run', 'main.nf', '-ansi-log', 'false',
                                     '-resume', '-with-trace', str(trace)], cwd=root,
                                    text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=120)
            assert result.returncode == 0, result.stdout
            with trace.open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            assert len(rows) == 1 and rows[0]['status'] == expected, rows
            hashes.append(rows[0]['hash'])
            print(f"{label}: {rows[0]['status']} {hashes[-1]}", flush=True)
        assert hashes[0] == hashes[1] and hashes[1] != hashes[2], hashes
    print('R-source cache invalidation passed.')


if __name__ == '__main__':
    main()
