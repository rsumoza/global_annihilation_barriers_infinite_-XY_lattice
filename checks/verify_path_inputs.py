#!/usr/bin/env python3
"""Exact input and local-geometry audit for the archived uniform XY upper path.

This checks actual binary payloads and rational consequences. It does not rerun
the old interval residual/Hessian/path/separator chain. Two geometry enclosures
are consumed from that authenticated antecedent, with that dependency explicit.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import struct
import sys
import zipfile


def need(condition, message):
    if not condition:
        raise ValueError(message)


def exact_json(payload):
    return json.loads(payload, parse_float=F)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def manifest(archive, prefix):
    rows = archive.read(prefix + 'FILES.sha256').decode().splitlines()
    return {line.split(maxsplit=1)[1].lstrip('*'): line.split(maxsplit=1)[0]
            for line in rows if line.strip()}


def checked_read(archive, prefix, relative, expected, records):
    payload = archive.read(prefix + relative)
    digest = hashlib.sha256(payload).hexdigest()
    need(expected.get(relative) == digest, 'Manifest mismatch: ' + relative)
    records[relative] = {'sha256': digest, 'bytes': len(payload)}
    return payload


def check_array(payload, side, nodes, label):
    need(sys.byteorder == 'little' and struct.calcsize('d') == 8,
         'Requires little-endian IEEE binary64 host')
    need(len(payload) == nodes * side * side * 8, 'Wrong array size: ' + label)
    view = memoryview(payload).cast('d')
    maximum = 0.0
    for value in view:
        need(math.isfinite(value), 'Nonfinite value: ' + label)
        maximum = max(maximum, abs(value))
    for k in range(nodes):
        offset = k * side * side
        for j in range(side):
            need(view[offset + j] == view[offset + (side - 1) * side + j] == 0,
                 'Nonzero row boundary: ' + label)
        for i in range(side):
            need(view[offset + i * side] == view[offset + i * side + side - 1] == 0,
                 'Nonzero column boundary: ' + label)
    return F.from_float(maximum)


def exact_gram(payload):
    need(len(payload) == 8192 * 8, 'Wrong two-mode payload')
    numbers = list(struct.unpack('<8192d', payload))
    need(all(math.isfinite(v) for v in numbers), 'Nonfinite reaction mode')
    q = [[F.from_float(v) for v in numbers[k * 4096:(k + 1) * 4096]] for k in range(2)]
    for k in range(2):
        for i in range(64):
            need(q[k][i] == q[k][63 * 64 + i] == q[k][i * 64] == q[k][i * 64 + 63] == 0,
                 'Reaction mode boundary not zero')
    g00 = g01 = g11 = F(0)
    for i in range(64):
        for j in range(64):
            v = i * 64 + j
            for t in ([v + 64] if i < 63 else []) + ([v + 1] if j < 63 else []):
                d0, d1 = q[0][t] - q[0][v], q[1][t] - q[1][v]
                g00 += d0 * d0
                g01 += d0 * d1
                g11 += d1 * d1
    det = g00 * g11 - g01 * g01
    need(det > 0, 'Reaction Gram matrix not positive definite')
    return [[g00, g01], [g01, g11]], [[g11 / det, -g01 / det], [-g01 / det, g00 / det]]


def interval_contains(interval, value):
    return F(interval[0]) <= value <= F(interval[1])


def ceil_fraction(value):
    return -(-value.numerator // value.denominator)


def audit(connection_path, optimality_path):
    connection_records, optimality_records = {}, {}
    cp, op = 'xy_global_connection/', 'xy_optimality_extension/'
    with zipfile.ZipFile(connection_path) as connection, zipfile.ZipFile(optimality_path) as optimality:
        cm, om = manifest(connection, cp), manifest(optimality, op)
        read_c = lambda rel: checked_read(connection, cp, rel, cm, connection_records)
        read_o = lambda rel: checked_read(optimality, op, rel, om, optimality_records)

        base_h = read_c('distribution/stage2/common/base_h.raw')
        h_min = check_array(base_h, 2048, 1, 'minimum corrector')
        del base_h
        saddle_h = read_c('distribution/stage2/saddle/inputs/h.raw')
        h_saddle = check_array(saddle_h, 4096, 1, 'saddle corrector')
        del saddle_h

        paths, maximum_z = [], F(0)
        for name, count in [('leg_return', 8), ('leg_annihilation', 19)]:
            cfg = list(map(int, read_c('distribution/connections/' + name + '.cfg').split()))
            need(cfg == [count, 256, 12], 'Wrong path header: ' + name)
            payload = read_c('distribution/connections/' + name + '.z.raw')
            zmax = check_array(payload, 256, count, name)
            if name == 'leg_return':
                need(all(v == 0 for v in memoryview(payload[:256 * 256 * 8]).cast('d')),
                     'Return leg does not start at the exact base reference')
            maximum_z = max(maximum_z, zmax)
            paths.append({'name': name, 'nodes': count, 'segments': count - 1,
                          'side': 256, 'max_abs_exact': zmax})
            del payload

        endpoint = exact_json(read_c('results_provided/connections/leg_annihilation.json'))
        need(endpoint['compact_integer_lift_max'] == 0 and endpoint['monotone_continuation_to_constant'],
             'Archived endpoint is not certified to contract to the zero lift')

        # Phi is in [-pi,pi] on this half-integer lattice. Core edges touching
        # |x|_infinity<=32 number <=4*64^2 and |Dphi|<=2*pi. Beyond that box,
        # |Dphi|<=4L/r^2=48/r^2; shell edge count <=64r. Sum r^{-3} by an
        # integral from 31 to infinity. These deliberately coarse bounds suffice.
        phi_square_upper = F(65536) * F(22, 7) ** 2 + F(73728, 31 ** 2)
        need(phi_square_upper < 1000 ** 2, 'Analytic reference norm bound failed')
        phi_norm = F(1000)
        rho_min, rho_saddle = F('0.017422'), F('0.012825')
        min_certificate = exact_json(read_c('results_provided/prerequisites/baseline.json'))
        saddle_certificate = exact_json(read_c('results_provided/prerequisites/saddle.json'))
        need(F(min_certificate['rho']) < rho_min and F(saddle_certificate['rho']) < rho_saddle,
             'Rational Newton-radius majorants fail')
        need(h_min < F(1, 5) and h_saddle < F(1, 2), 'Candidate phase bounds fail')
        min_norm = phi_norm + 3 * 2048 * h_min + rho_min
        saddle_norm = phi_norm + 3 * 4096 * h_saddle + rho_saddle
        explicit_vertex_norm = phi_norm + 3 * 2048 * h_min + 3 * 256 * maximum_z
        input_radius = F(1, 50)
        input_norm = phi_norm + 3 * 2048 * h_min + input_radius
        uniform_integer_norm = ceil_fraction(max(min_norm, saddle_norm, explicit_vertex_norm, input_norm)) + 1

        # The half-integer coordinates of the two correctors lie strictly inside
        # radius 4096. The extra two lattice units below cover both a unit
        # edge and the half-lattice translation when integer box coordinates
        # are used instead of the physical half-integer coordinates.
        r0, alpha = 4096, F('0.018')
        dipole_increment = F(12, (r0 - 12 - 2) ** 2)
        need(dipole_increment + max(rho_min, rho_saddle) < alpha,
             'Stationary exterior increment bound fails')
        ellipticity_lower = F(9999, 10000)
        need(1 - alpha ** 2 / 6 > ellipticity_lower, 'Exterior ellipticity bound fails')

        q_payload = read_c('distribution/stage2/barrier/inputs/q.raw')
        gram, gram_inverse = exact_gram(q_payload)
        caps = exact_json(read_o('results_provided/caps.json'))
        need(all(interval_contains(caps['G'][i][j], gram[i][j]) for i in range(2) for j in range(2)),
             'Reaction mode Gram differs from the capped separator')
        need(gram_inverse[0][0] < 4 and gram_inverse[1][1] < 4,
             'Coordinate-functional norm bound fails')
        geometry = exact_json(read_o('antecedent_results/separator/geometry.json'))
        center = next(c for c in geometry['centers'] if c['grid'] == [0, 0])
        need(center['id'] == 54, 'Central chart id differs')
        central_chart = exact_json(read_o('antecedent_results/separator/charts/c054/chart.json'))
        need(F(central_chart['r']) >= F('0.08') and F(central_chart['R']) >= F('0.25'),
             'Central chart too small')
        center_error = F(1, 10 ** 12)
        need(all(max(abs(F(v)) for v in axis) < center_error for axis in center['a']),
             'Central coordinate error too large')
        need(F(center['Y_norm_squared'][1]) < center_error ** 2, 'Central transverse error too large')
        coordinate_bound = 2 * input_radius
        need(coordinate_bound + center_error < F(central_chart['r']), 'Coordinate chart membership fails')
        need(input_radius + center_error < F(central_chart['R']), 'Transverse chart membership fails')
        disk_radius = F('0.6599')
        need(2 * coordinate_bound < disk_radius, 'Disk membership fails')
        cut_margins = []
        for cap in caps['caps']:
            normal_l1 = sum((max(abs(F(v)) for v in axis) for axis in cap['normal']), F(0))
            margin = F(cap['offset'][0]) - normal_l1 * coordinate_bound
            need(margin > 0, 'Cut membership fails: ' + cap['name'])
            cut_margins.append({'name': cap['name'], 'margin_lower': margin})
        energy_gap = F('0.00567867')
        segment_energy_upper = (input_radius + rho_min) ** 2 / 2
        need(segment_energy_upper < energy_gap, 'Input-to-minimum segment might reach saddle height')

        return encode({
            'status': 'UPPER_PATH_INPUT_AND_LOCAL_CLASS_AUDIT_PASS',
            'full_interval_chain_reexecuted': False,
            'antecedent_dependency': 'Authenticated minimum, saddle, path, and capped-separator certificates remain prerequisites; this audit checks payloads and rational consequences.',
            'conditional_on_antecedent_geometry_enclosures': True,
            'model': 'Homogeneous unit-coupling exact cosine on the square lattice',
            'physical_coordinates': 'x=i-(n-1)/2+6, y=j-(n-1)/2',
            'reference': 'phi12(x,y)=Arg(x+i*y)-Arg(x-12+i*y)',
            'minimum_corrector_max_exact': h_min, 'saddle_corrector_max_exact': h_saddle,
            'compact_paths': paths,
            'upper_path_sections': ['m_star to m_hat: convex straight segment',
                'leg_return: 7 straight segments', 'c_plus to s_star: Hermite straight segment',
                's_star to c_minus: Hermite straight segment', 'leg_annihilation: 18 straight segments',
                'endpoint to zero real lift: straight monotone contraction'],
            'implicit_vertices': ['m_star', 's_star'],
            'reference_norm_square_upper': phi_square_upper, 'reference_norm_upper': phi_norm,
            'compact_norm_bound': '||Dh||_2 <= sqrt(8)||h||_2 <= 3*n*max_abs(h)',
            'minimum_exact_norm_upper': min_norm, 'saddle_exact_norm_upper': saddle_norm,
            'explicit_path_vertex_norm_upper': explicit_vertex_norm,
            'input_class_norm_upper': input_norm, 'uniform_path_and_critical_norm_upper_integer': uniform_integer_norm,
            'candidate_core_safe_radius': r0, 'candidate_exterior_increment_upper': dipole_increment,
            'exact_stationary_exterior_increment_upper': alpha,
            'stationary_exterior_ellipticity_lower': ellipticity_lower,
            'ellipticity_justification': 'sin(t)/t >= 1-t^2/6 for |t|<=0.018',
            'input_class': {'description': 'Open X-ball B_X(m_hat,1/50)', 'radius': input_radius,
                'contains_exact_minimum': rho_min < input_radius, 'central_chart_id': center['id'],
                'coordinate_functional_norm_upper': 2, 'coordinate_absolute_bound': coordinate_bound,
                'central_coordinate_error_upper': center_error, 'central_transverse_error_upper': center_error,
                'inside_chart_disk_and_both_cuts': True, 'cut_margins_lower': cut_margins,
                'energy_above_exact_minimum_upper_on_segment': segment_energy_upper,
                'certified_saddle_gap_lower_used': energy_gap, 'strict_energy_margin': energy_gap - segment_energy_upper},
            'reaction_gram_exact': gram, 'reaction_gram_inverse_exact': gram_inverse,
            'input_sha256': {'connection': connection_records, 'optimality': optimality_records}
        })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--connection-zip', type=Path, required=True)
    parser.add_argument('--optimality-zip', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.connection_zip, args.optimality_zip)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'],
                      'M': result['uniform_path_and_critical_norm_upper_integer'],
                      'input_class': result['input_class']['description'],
                      'output': str(args.output)}, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError, zipfile.BadZipFile) as error:
        print('REJECT: ' + str(error), file=sys.stderr)
        sys.exit(2)
