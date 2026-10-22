#!/usr/bin/env python3
"""Two intentionally invalid geometry ledgers must be rejected, also under -O."""
from argparse import ArgumentParser
from pathlib import Path
import json
import subprocess
import sys
import tempfile


def main():
    ap = ArgumentParser()
    ap.add_argument('--package', type=Path, required=True)
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    checks = []
    with tempfile.TemporaryDirectory(prefix='geometry_negative_') as name:
        directory = Path(name)
        for kind in ('excluded_coordinate', 'missing_neighbor'):
            data = json.loads(args.geometry.read_text())
            if kind == 'excluded_coordinate':
                data['centers'][0]['a'][0] = ['0', '0']
            else:
                data['neighbors'].pop()
            mutated = directory / (kind + '.json')
            mutated.write_text(json.dumps(data))
            output = directory / (kind + '_report.json')
            command = [sys.executable, '-O', str(Path(__file__).with_name('check_geometry_exact.py')),
                       '--package', str(args.package.resolve()), '--geometry', str(mutated), '--output', str(output)]
            result = subprocess.run(command, text=True, capture_output=True)
            if result.returncode == 0 or output.exists():
                raise RuntimeError('negative control falsely accepted: ' + kind)
            expected = {'excluded_coordinate': 'exact value outside interval: a.0.0',
                        'missing_neighbor': 'full neighbor ledger differs'}[kind]
            if expected not in result.stderr:
                raise RuntimeError('control rejected for unexpected reason: ' + result.stderr)
            checks.append({'mutation': kind, 'returncode': result.returncode,
                           'pass_output_written': output.exists(),
                           'rejection': result.stderr.splitlines()[-1]})
    report = {'status': 'EXACT_GEOMETRY_NEGATIVE_CONTROLS_PASS', 'python_optimized': True,
              'controls': checks, 'scope': 'Two intentional output mutations, not exhaustive parser fuzzing.'}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
