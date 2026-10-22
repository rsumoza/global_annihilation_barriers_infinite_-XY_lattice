#!/usr/bin/env python3
"""Verify unchanged original archives before extraction; never executes them."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets-dir', type=Path, required=True)
    parser.add_argument('--manifest', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'release-assets.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    results = []
    for item in manifest['assets']:
        name = item['name']
        if Path(name).name != name:
            raise ValueError('Asset manifest contains a nonlocal filename')
        path = args.assets_dir / name
        if not path.is_file():
            raise ValueError('Missing required original: ' + name)
        size = path.stat().st_size
        if size != item['bytes']:
            raise ValueError('Byte length mismatch: ' + name)
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(chunk)
        value = digest.hexdigest()
        if value != item['sha256']:
            raise ValueError('SHA-256 mismatch: ' + name)
        results.append({'item': item['item'], 'name': name,
                        'bytes': size, 'sha256': value})
    report = {'status': 'PASS', 'scope': 'Original archive identity only; no extraction or mathematical replay.',
              'assets': results}
    text = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
