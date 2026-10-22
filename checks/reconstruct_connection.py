#!/usr/bin/env python3
"""Reconstruct the exact frozen connection ZIP from the five supplied volumes.

Does not modify the volumes or extract/run any contained checking source.
The completed file is installed atomically only after its full hash matches.
"""
from argparse import ArgumentParser
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import os
import tempfile

EXPECTED_BYTES = 441730856
EXPECTED_SHA256 = '3a4163d8b24c21f95acbdf502ceaed46e967180c21abb363d3436ae5e2afb042'
VOLUME_HASHES = [
    'db52c9b733276c6e2ed4b2f1c4af3e125196239337b06734dda0c05f2b322c9b',
    '3c28659c4d13cc664e3b0e4f2970eaa84231450f5b4e8134b895c9717876a8f9',
    '201ce9a0c73d4ccafda83ea551f60bfdc9fefdd716803926b0072491a831660c',
    '9842cfca8694db2d762c21edb4ddd3fc3b051eaf8c6fe64a945a18c6eb4b8b06',
    '0e03784aab5da0aa63268121225c3ad5713fdcb954fe508fd23df64bc0bddf60',
]

def need(value, reason):
    if not value:
        raise RuntimeError(reason)

def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--volumes', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    records = []
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=args.output.parent, delete=False) as target:
            temporary = Path(target.name)
            for number, expected in enumerate(VOLUME_HASHES, 1):
                path = args.volumes / f'S3_Conexion_{number:02}_de_05.zip'
                need(sha(path) == expected, 'Volume identity mismatch: ' + path.name)
                member = f'payload/S3_{number:03}.part'
                digest = hashlib.sha256()
                size = 0
                with ZipFile(path) as archive:
                    need(sorted(archive.namelist()) == sorted(['README_VOLUME.txt', member]),
                         'Volume member ledger differs: ' + path.name)
                    with archive.open(member) as stream:
                        while block := stream.read(4 * 1024 * 1024):
                            target.write(block)
                            digest.update(block)
                            size += len(block)
                records.append({'volume': path.name, 'sha256': expected,
                                'fragment_bytes': size, 'fragment_sha256': digest.hexdigest()})
            target.flush()
            os.fsync(target.fileno())
        need(temporary.stat().st_size == EXPECTED_BYTES and sha(temporary) == EXPECTED_SHA256,
             'Reconstructed archive does not match the frozen original')
        os.replace(temporary, args.output)
        temporary = None
        need(args.output.stat().st_size == EXPECTED_BYTES and sha(args.output) == EXPECTED_SHA256,
             'Installed archive identity differs')
        print(json.dumps({'status': 'EXACT_CONNECTION_ARCHIVE_RECONSTRUCTED',
                          'bytes': EXPECTED_BYTES, 'sha256': EXPECTED_SHA256,
                          'output': str(args.output), 'volumes': records}, indent=2))
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)

if __name__ == '__main__':
    main()
