#!/usr/bin/env python3
"""Compare freshly rebuilt separator outputs with the frozen record.

This is a regression/binding check, not a second arithmetic implementation.
Only path/timing dependent fields are removed, explicitly below.
"""
from pathlib import Path
from fractions import Fraction
from argparse import ArgumentParser
import hashlib
import json

ROOT = Path(__file__).resolve().parent
parser = ArgumentParser(description=__doc__)
parser.add_argument('--package', type=Path, default=ROOT / 'xy_global_separator')
parser.add_argument('--fresh', type=Path, default=ROOT / 'fresh_results')
parser.add_argument('--output', type=Path, default=ROOT / 'fresh_comparison.json')
args = parser.parse_args()
OLD = args.package.resolve() / 'results_provided'
NEW = args.fresh.resolve()

def require(value, reason):
    if not value:
        raise RuntimeError(reason)

def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()

def read(path):
    return json.loads(path.read_text(), parse_float=Fraction)

checked = []
for name in ('baseline.json', 'geometry.json', 'grid.json'):
    require((OLD / name).read_bytes() == (NEW / name).read_bytes(), 'changed fixed output: ' + name)
    checked.append(name)
for directory in sorted((OLD / 'charts').iterdir()):
    require(directory.is_dir(), 'unexpected chart entry')
    for name in ('chart.json', 'coefficients.txt', 'polynomials.json'):
        relative = Path('charts') / directory.name / name
        require((OLD / relative).read_bytes() == (NEW / relative).read_bytes(), 'changed chart: ' + str(relative))
        checked.append(str(relative))
require(len(checked) == 330, 'chart file ledger is incomplete')

old_boundary = read(OLD / 'boundary.json')
new_boundary = read(NEW / 'boundary.json')
for obj in (old_boundary, new_boundary):
    obj.pop('elapsed_seconds', None)
    obj.pop('input_hashes', None)
    for chart in obj['charts']:
        chart.pop('path', None)
require(old_boundary == new_boundary, 'boundary numerical/topological content differs')
old_separation = read(OLD / 'separation.json')
new_separation = read(NEW / 'separation.json')
for obj in (old_separation, new_separation):
    obj.pop('source_sha256', None)
require(old_separation == new_separation, 'separation content differs')

certificate = read(NEW / 'certificate.json')
require(certificate['status'] == 'FRESH_GLOBAL_SEPARATOR_PASS', 'fresh supervisor did not pass')
require(certificate['manifest_sha256'] == '991c744eab7e812d629662a0c0bc5cbc177e9dd2e9daa67caeefab7a814f941b', 'scientific pin changed')
require(certificate['supervisor_sha256'] == 'd4646d988e538075b51d0c0bdd2d970976d1ccb53156ad52d414ec87f5c0095c', 'supervisor pin changed')
for name, expected in certificate['output_sha256'].items():
    require(sha(NEW / name) == expected, 'fresh output hash mismatch: ' + name)

result = {
    'status': 'FRESH_SEPARATOR_REPRODUCTION_COMPARISON_PASS',
    'scope': 'Fresh unchanged frozen software; byte/numerical comparison only, not an independent interval backend',
    'byte_identical_scientific_outputs': checked,
    'byte_identical_output_count': len(checked),
    'boundary_equal_after_explicit_path_timing_normalization': True,
    'separation_equal_after_explicit_source_hash_normalization': True,
    'fresh_output_hashes_checked': len(certificate['output_sha256']),
    'certificate_sha256': sha(NEW / 'certificate.json'),
    'original_sources_changed': False,
    'connection_reexecuted': False,
    'optimality_reexecuted': False,
}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
