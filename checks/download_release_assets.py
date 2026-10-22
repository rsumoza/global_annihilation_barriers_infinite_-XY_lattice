#!/usr/bin/env python3
"""Download original S2-S4 Release assets; verify exact bytes before installing.

The configured destinations are not evidence that the Release is published.
Uses Python's standard library. Does not extract or execute downloaded files.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


def identity(path):
    h = hashlib.sha256()
    count = 0
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            count += len(chunk)
            h.update(chunk)
    return count, h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets-dir', type=Path)
    parser.add_argument('--manifest', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'release-assets.json')
    parser.add_argument('--list', action='store_true', help='List URLs without network access')
    args = parser.parse_args()
    if not args.list and args.assets_dir is None:
        parser.error('--assets-dir is required unless --list is used')
    manifest = json.loads(args.manifest.read_text())
    root = manifest['repository_url'].rstrip('/')
    if not root.startswith('https://github.com/') or root.endswith('.git'):
        raise ValueError('Expected an HTTPS GitHub repository URL without .git')
    tag = manifest['planned_tag']
    if not args.list:
        args.assets_dir.mkdir(parents=True, exist_ok=True)
    for item in manifest['assets']:
        name = item['name']
        if not name or Path(name).name != name or '\\' in name:
            raise ValueError('Manifest contains an invalid local filename')
        url = root + '/releases/download/' + quote(tag, safe='') + '/' + quote(name, safe='')
        if args.list:
            print(item['item'], item['bytes'], item['sha256'], url)
            continue
        target = args.assets_dir / name
        expected = (item['bytes'], item['sha256'])
        if target.exists():
            if identity(target) != expected:
                raise ValueError('Existing file differs; preserve or move it before retrying: ' + str(target))
            print(item['item'] + ': existing verified original ' + name)
            continue
        temporary = None
        try:
            request = Request(url, headers={'User-Agent': 'XY-SIADS-reproducibility-v7'})
            with urlopen(request, timeout=60) as response:
                with tempfile.NamedTemporaryFile(dir=args.assets_dir, prefix='.' + name,
                                                 suffix='.part', delete=False) as out:
                    temporary = Path(out.name)
                    h = hashlib.sha256()
                    count = 0
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        count += len(chunk)
                        if count > item['bytes']:
                            raise ValueError('Response exceeds expected size: ' + name)
                        out.write(chunk)
                        h.update(chunk)
                    out.flush()
                    os.fsync(out.fileno())
            if (count, h.hexdigest()) != expected:
                raise ValueError('Byte count or SHA-256 mismatch: ' + name)
            if target.exists():
                raise FileExistsError('Destination appeared during download: ' + str(target))
            os.link(temporary, target)
            print(item['item'] + ': downloaded and verified ' + name)
        except (HTTPError, URLError) as exc:
            raise SystemExit('Release asset unavailable. Verify that the planned Release is published, '
                             'publicly accessible, and contains the canonical asset names. ' + str(exc)) from exc
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
