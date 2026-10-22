#!/usr/bin/env python3
"""Independent all-entry exact-dyadic audit of separator geometry.

Uses only Python's standard library. Does not import the supplied engine or
its checkers. Input IEEE binary64 numbers are exact rationals. Dirichlet inner
products are calculated with the full-lattice degree-four Laplacian applied
to zero extensions, rather than reproducing C++ interval edge accumulation.

Scope: all reported geometry enclosures and scalar maxima, not interval
analytic bounds for the infinite minimum or the chart energies.
"""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import struct
import time

PIN = '991c744eab7e812d629662a0c0bc5cbc177e9dd2e9daa67caeefab7a814f941b'


def require(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read_dyadic(path, count):
    data = path.read_bytes()
    require(len(data) == 8 * count, 'binary64 array length: ' + str(path))
    values = struct.unpack('<' + str(count) + 'd', data)
    require(all(math.isfinite(v) for v in values), 'nonfinite input: ' + str(path))
    return [v.as_integer_ratio() for v in values]


def laplace_zero(values, n):
    """Full-lattice degree-four graph Laplacian of the zero extension."""
    result = [4 * v for v in values]
    for i in range(n):
        for j in range(n):
            k = i * n + j
            if i:
                result[k] -= values[k - n]
            if i + 1 < n:
                result[k] -= values[k + n]
            if j:
                result[k] -= values[k - 1]
            if j + 1 < n:
                result[k] -= values[k + 1]
    return result


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def check_interval(value, interval, name, checks):
    require(isinstance(interval, list) and len(interval) == 2, 'interval shape: ' + name)
    low, high = map(F, interval)
    require(low <= value <= high, 'exact value outside interval: ' + name)
    checks.append(name)


def check_sqrt_interval(square, interval, name, checks):
    low, high = map(F, interval)
    require(square >= 0 and 0 <= low <= high, 'invalid square/distance: ' + name)
    require(low * low <= square <= high * high, 'exact distance outside interval: ' + name)
    checks.append(name)


def run(package, geometry_path):
    start = time.perf_counter()
    src = package / 'distribution'
    manifest = src / 'MANIFEST.sha256'
    require(sha(manifest) == PIN, 'scientific manifest pin mismatch')
    bindings = {}
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        bindings[name] = digest
    used = ['common/q.raw', 'geometry/cfg', 'grid.json']
    rows = (src / 'geometry/cfg').read_text().splitlines()
    n, count = map(int, rows[0].split())
    require((n, count) == (128, 109), 'unexpected geometry dimensions')
    records = []
    for expected, row in enumerate(rows[1:]):
        token = row.split()
        require(len(token) == 4 and int(token[0]) == expected, 'configuration ledger')
        name = token[3]
        require(Path(name).name == name and '..' not in name, 'unsafe center name')
        records.append({'id': expected, 'grid': [int(token[1]), int(token[2])], 'file': name})
        used.append('centers/' + name)
    require(len(records) == count, 'configuration record count')
    for name in used:
        require(name in bindings and sha(src / name) == bindings[name], 'bound input mismatch: ' + name)
    grid = json.loads((src / 'grid.json').read_text())
    require(grid['n'] == n and len(grid['centers']) == count, 'grid dimensions')
    for r, g in zip(records, grid['centers']):
        require((r['id'], r['grid'], r['file']) == (g['id'], g['grid_index'], g['file']), 'grid/cfg disagreement')

    rawq = read_dyadic(src / 'common/q.raw', 8192)
    q_small = [rawq[:4096], rawq[4096:]]
    z_small = [read_dyadic(src / 'centers' / r['file'], n * n) for r in records]
    # Every denominator is a power of two. A single exact integer unit supports
    # all arrays; the largest denominator is found from the actual input bits.
    denominator = max(d for num, d in rawq)
    for values in z_small:
        denominator = max(denominator, max(d for num, d in values))
    require(denominator & (denominator - 1) == 0, 'non-dyadic denominator')
    scale2 = denominator * denominator
    q = []
    offset = (n - 64) // 2
    for v in q_small:
        full = [0] * (n * n)
        for i in range(64):
            for j in range(64):
                num, den = v[i * 64 + j]
                require(i not in (0, 63) and j not in (0, 63) or num == 0, 'reaction mode nonzero boundary')
                full[(i + offset) * n + j + offset] = num * (denominator // den)
        q.append(full)
    z = [[num * (denominator // den) for num, den in values] for values in z_small]
    del rawq, q_small, z_small
    for values in z:
        require(all(values[k] == 0 for k in range(n)) and
                all(values[(n - 1) * n + k] == 0 for k in range(n)) and
                all(values[k * n] == values[k * n + n - 1] == 0 for k in range(n)),
                'center nonzero boundary; supplied inner-domain formula would omit exterior bonds')
    lq = [laplace_zero(values, n) for values in q]
    G = [[F(dot(q[k], lq[l]), scale2) for l in range(2)] for k in range(2)]
    require(G[0][1] == G[1][0], 'Gram symmetry')
    det = G[0][0] * G[1][1] - G[0][1] * G[1][0]
    require(G[0][0] > 0 and det > 0, 'non-positive exact Gram matrix')
    Gi = [[G[1][1] / det, -G[0][1] / det], [-G[1][0] / det, G[0][0] / det]]
    target = json.loads(geometry_path.read_text())
    require(target['status'] == 'ATLAS_PROJECTION_AND_NEIGHBOR_GEOMETRY_PASS', 'geometry status')
    require((target['n'], target['count']) == (n, count), 'output dimensions')
    checks = []
    for k in range(2):
        for l in range(2):
            check_interval(G[k][l], target['G'][2 * k + l], f'G.{k}.{l}', checks)
            check_interval(Gi[k][l], target['Gi'][2 * k + l], f'Gi.{k}.{l}', checks)
    g = [[F(dot(values, lq[k]), scale2) for k in range(2)] for values in z]
    lz = [laplace_zero(values, n) for values in z]
    norm = [F(dot(values, lap), scale2) for values, lap in zip(z, lz)]
    a = [[sum(Gi[k][l] * gl[l] for l in range(2)) for k in range(2)] for gl in g]
    y2 = [norm[k] - sum(g[k][l] * a[k][l] for l in range(2)) for k in range(count)]
    centers = target['centers']
    require(len(centers) == count, 'output center count')
    max_y_upper = F(target['maximum_center_Y_norm_upper'])
    max_coord_upper = F(target['maximum_nominal_coordinate_error_upper'])
    require(max_y_upper >= 0 and max_coord_upper >= 0, 'negative scalar upper bound')
    max_exact_coordinate_error = F(0)
    for k, (rec, center) in enumerate(zip(records, centers)):
        require(rec == {key: center[key] for key in ('id', 'grid', 'file')}, 'center output identity')
        for l in range(2):
            check_interval(a[k][l], center['a'][l], f'a.{k}.{l}', checks)
            error = abs(a[k][l] - F(3 * rec['grid'][l], 25))
            max_exact_coordinate_error = max(max_exact_coordinate_error, error)
            require(error <= max_coord_upper, 'coordinate maximum enclosure')
        require(y2[k] >= 0, 'negative exact Y norm squared')
        check_interval(y2[k], center['Y_norm_squared'], f'Y2.{k}', checks)
        require(y2[k] <= max_y_upper * max_y_upper, 'Y norm maximum enclosure')
    pairs = [(i, j) for i in range(count) for j in range(i + 1, count)
             if max(abs(records[i]['grid'][k] - records[j]['grid'][k]) for k in range(2)) <= 1]
    require(len(pairs) == 380 and target['neighbor_pairs'] == len(pairs), 'neighbor count')
    require([(v['i'], v['j']) for v in target['neighbors']] == pairs, 'full neighbor ledger differs')
    max_distance_upper = F(target['maximum_neighbor_Y_distance_upper'])
    require(max_distance_upper >= 0, 'negative neighbor maximum')
    max_square = F(0)
    for (i, j), row in zip(pairs, target['neighbors']):
        norm_difference = norm[i] + norm[j] - F(2 * dot(z[i], lz[j]), scale2)
        gd = [g[i][k] - g[j][k] for k in range(2)]
        projected_square = norm_difference - sum(gd[k] * Gi[k][l] * gd[l] for k in range(2) for l in range(2))
        check_sqrt_interval(projected_square, row['distance'], f'distance.{i}.{j}', checks)
        require(projected_square <= max_distance_upper * max_distance_upper, 'neighbor maximum enclosure')
        max_square = max(max_square, projected_square)
    require(len(checks) == 8 + 3 * count + len(pairs), 'incomplete enclosure count')
    return {
        'status': 'INDEPENDENT_ALL_GEOMETRY_EXACT_DYADIC_PASS',
        'scope': 'All fixed geometry entries and reported maxima; excludes minimum/chart analytic estimates and global topology.',
        'source_package': str(package.resolve()),
        'geometry_file': str(geometry_path.resolve()),
        'geometry_sha256': sha(geometry_path),
        'scientific_manifest_sha256': PIN,
        'input_sha256': {name: bindings[name] for name in used},
        'checker_sha256': sha(Path(__file__)),
        'method': 'Exact binary64 dyadic decode; arbitrary-size integer degree-four zero-extension Laplacian; rational Gram inverse and squared-distance comparisons.',
        'common_dyadic_denominator_log2': denominator.bit_length() - 1,
        'raw_arrays_checked': count + 2,
        'reaction_mode_boundary_checks': True,
        'center_boundary_checks': True,
        'includes_zero_extension_exterior_degree': True,
        'gram_and_inverse_enclosures': 8,
        'center_coordinate_enclosures': 2 * count,
        'center_Y_squared_enclosures': count,
        'neighbor_distance_enclosures': len(pairs),
        'total_scalar_interval_enclosures': len(checks),
        'scalar_maxima_verified': ['maximum_neighbor_Y_distance_upper', 'maximum_center_Y_norm_upper', 'maximum_nominal_coordinate_error_upper'],
        'maximum_exact_neighbor_squared': str(max_square),
        'maximum_exact_coordinate_error': str(max_exact_coordinate_error),
        'duration_seconds': time.perf_counter() - start,
    }


def main():
    ap = ArgumentParser()
    ap.add_argument('--package', type=Path, required=True)
    ap.add_argument('--geometry', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    result = run(args.package, args.geometry)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('status', 'total_scalar_interval_enclosures', 'duration_seconds')}, indent=2))


if __name__ == '__main__':
    main()
