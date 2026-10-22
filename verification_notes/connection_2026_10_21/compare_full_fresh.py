#!/usr/bin/env python3
"""Validate the fresh global record and compare its scientific outputs.

All paths/timings removed from comparisons are named explicitly. This
regression check is not a second interval implementation or an analytic proof.
"""
from argparse import ArgumentParser
from pathlib import Path
from fractions import Fraction
import hashlib
import json

HERE = Path(__file__).resolve().parent
PINS = {
    'supervisor_sha256': '564451be670e23dd08d382f0f71fc9bf238a39fb7003db916c3f953720e840c3',
    'extension_manifest_sha256': '248f0614524832610b07fd21b5d6b9eaf670d1846d6032aa22a51755a5c95972',
    'separator_manifest_sha256': '991c744eab7e812d629662a0c0bc5cbc177e9dd2e9daa67caeefab7a814f941b',
    'connection_manifest_sha256': 'ba100e9b7b8763bcd9ce1bcb7ee1a50a1f76669d2be92966ecbe6f6208feb27c',
}
CHANGED = {42, 43, 44, 53, 55}

def need(value, reason):
    if not value:
        raise RuntimeError(reason)

def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def read(path):
    return json.loads(path.read_text(), parse_float=Fraction)

def normalize_boundary(value):
    value.pop('elapsed_seconds', None)
    value.pop('input_hashes', None)
    for chart in value['charts']:
        chart.pop('path', None)
    return value

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--fresh', type=Path, default=HERE / 'full_fresh')
    parser.add_argument('--connection', type=Path, required=True)
    parser.add_argument('--optimality', type=Path, required=True)
    parser.add_argument('--historical-certificate', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE / 'full_comparison.json')
    args = parser.parse_args()
    new = args.fresh.resolve()
    old = args.optimality.resolve() / 'results_provided'
    c_old = args.connection.resolve() / 'results_provided'
    cert = read(new / 'certificate.json')
    need(cert['status'] == 'FRESH_OPTIMAL_ANNIHILATION_BARRIER_PASS', 'fresh global result absent')
    need(cert['fresh_results_only'] and cert['fresh_extension_computations'] and
         cert['antecedent_execution_mode'] == 'fresh_subprocesses', 'resumed execution is not a fresh run')
    for key, value in PINS.items():
        need(cert[key] == value, 'frozen identity changed: ' + key)
    historical = read(args.historical_certificate)
    mathematical_keys = [
        'scope', 'energy_normalization', 'identity', 'saddle_energy_shared_by_exact_reflection',
        'global_optimality_proved', 'heteroclinic_orbit_proved', 'universal_XY_metastability_proved',
        'all_boundary_parts_verified', 'symmetry_restriction_on_perturbations',
        'barrier_lower_exact', 'barrier_upper_exact', 'reported_strict_bracket', 'disk_radius',
        'charts', 'coordinate_faces', 'cap_coercivity_safe', 'cap_radius', 'flags',
        *PINS,
    ]
    for key in mathematical_keys:
        need(cert[key] == historical[key], 'mathematical certificate changed: ' + key)
    checked_hashes = {}
    for path, expected in cert['output_sha256'].items():
        relative = Path(path)
        need(not relative.is_absolute() and '..' not in relative.parts, 'unsafe recorded output')
        need(sha(new / relative) == expected, 'fresh output hash mismatch: ' + path)
        checked_hashes[path] = expected
    nested = [
        ('separator', 'FRESH_GLOBAL_SEPARATOR_PASS'),
        ('connection', 'FRESH_GLOBAL_VARIATIONAL_SADDLE_CONNECTION_PASS'),
    ]
    nested_checks = {}
    for name, status in nested:
        folder = new / ('fresh_' + name)
        need(sha(folder / 'certificate.json') == cert['antecedent_certificate_sha256'][name],
             'antecedent certificate binding differs: ' + name)
        stage = read(folder / 'certificate.json')
        need(stage['status'] == status, 'antecedent rejected: ' + name)
        for path, expected in stage['output_sha256'].items():
            relative = Path(path)
            need(not relative.is_absolute() and '..' not in relative.parts, 'unsafe nested path')
            need(sha(folder / relative) == expected, 'nested hash mismatch: ' + name + '/' + path)
        nested_checks[name] = len(stage['output_sha256'])
    byte_identical = []
    connection_files = [
        'prerequisites/baseline.json', 'prerequisites/saddle.json',
        'prerequisites/barrier.json', 'prerequisites/rational.json',
        'connections/leg_return.json', 'connections/leg_annihilation.json',
        'connections/full_annihilation.json', 'connections/join_points.json',
        'global_rational.json',
    ]
    for path in connection_files:
        need((new / 'fresh_connection' / path).read_bytes() == (c_old / path).read_bytes(),
             'connection scientific bytes differ: ' + path)
        byte_identical.append('fresh_connection/' + path)
    for path in ['projection_plus.json', 'projection_minus.json',
                 'adjoint_plus.json', 'adjoint_minus.json',
                 'center_distances_plus.json', 'center_distances_minus.json',
                 'caps.json', 'domain.json', 'membership.json', 'grid.json',
                 'antecedent_input_bindings.json']:
        need((new / path).read_bytes() == (old / path).read_bytes(),
             'extension scientific bytes differ: ' + path)
        byte_identical.append(path)
    manifest = read(new / 'atlas_manifest.json')['charts']
    need([row['id'] for row in manifest] == list(range(109)), 'full atlas ledger changed')
    for row in manifest:
        cid = row['id']
        target = Path(row['path'])
        if not target.is_absolute():
            target = new / target
        expected = (old / 'enlarged_charts' if cid in CHANGED else
                    args.optimality / 'antecedent_results/separator/charts') / f'c{cid:03}' / 'polynomials.json'
        need(target.read_bytes() == expected.read_bytes(), f'atlas polynomial bytes differ: {cid}')
        byte_identical.append(f'atlas/c{cid:03}/polynomials.json')
    for path in ('boundary.json', 'boundary_old_boundary.json'):
        need(normalize_boundary(read(new / path)) == normalize_boundary(read(old / path)),
             'boundary numerical/geometric content differs: ' + path)
    result = {
        'status': 'FULL_FRESH_GLOBAL_REPRODUCTION_COMPARISON_PASS',
        'scope': 'Fresh unchanged frozen checking software plus explicit data/regression checks; not an independent interval engine',
        'resumed_results_used': False,
        'fresh_certificate_sha256': sha(new / 'certificate.json'),
        'fresh_output_hashes_checked': len(checked_hashes),
        'nested_output_hashes_checked': nested_checks,
        'historical_mathematical_fields_equal': mathematical_keys,
        'byte_identical_scientific_outputs': byte_identical,
        'byte_identical_output_count': len(byte_identical),
        'both_boundary_records_equal_after_path_timing_normalization': True,
        'original_frozen_sources_changed': False,
        'heteroclinic_orbit_certified': False,
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('byte_identical_scientific_outputs',
                                     'historical_mathematical_fields_equal')}, indent=2))

if __name__ == '__main__':
    main()
