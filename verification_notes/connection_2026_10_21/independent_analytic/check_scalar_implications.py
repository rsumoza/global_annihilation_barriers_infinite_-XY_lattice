#!/usr/bin/env python3
"""Independent exact scalar audit of the connection's analytic implications.

This does not import the supplied checkers or reverify the C++ interval bounds.
It recomputes Newton, energy, local-box and Hermite implications from the
specified local outputs, with exact Fraction parsing and explicit padding.
"""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
DEFAULT = HERE.parent / 'xy_global_connection/results_provided'
PAD = F('1e-18')


def need(value, reason):
    if not value:
        raise RuntimeError(reason)


def read(path):
    return json.loads(path.read_text(), parse_float=F)


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def audit(root):
    pre = root / 'prerequisites'
    con = root / 'connections'
    b = read(pre / 'baseline.json')
    s = read(pre / 'saddle.json')
    joins = read(con / 'join_points.json')
    ret = read(con / 'leg_return.json')
    ann = read(con / 'leg_annihilation.json')
    need(b['status'] == 'INDEPENDENT_CHAIN_PASS', 'minimum prerequisite')
    need(s['status'] == 'INFINITE_SADDLE_INDEX_ONE_PASS', 'saddle prerequisite')
    need(s['morse_index'] == 1 and not s['symmetry_restriction_on_perturbations'], 'full-space index prerequisite')
    need(ret['status'] == 'COMPACT_PATH_LEG_PASS_RELATIVE_REFERENCE', 'return leg')
    need(ann['status'] == 'COMPACT_TERMINAL_PATH_PASS_RELATIVE_REFERENCE', 'annihilation leg')
    d = F(s['delta']) + PAD
    rho = F(s['rho']) + PAD
    a = F(s['negative_direction_bound']) - PAD
    m = F(s['transverse_coercivity']) - PAD
    gamma = min(a, m)
    invariant_slack = gamma * rho - d - rho * rho / 2
    need(invariant_slack >= 0 and rho < gamma, 'frozen Newton ball')
    astar, mstar = a - rho, m - rho
    need(min(astar, mstar) > F('0.057'), 'exact saddle gap')
    need(rho < F(s['candidate_branch_margin']) - PAD, 'charges along correction')
    eps = min(rho * rho / 2, d * rho / 2 + rho**3 / 12)
    slow = F(s['candidate_relative_to_reference'][0]) - PAD - eps
    shigh = F(s['candidate_relative_to_reference'][1]) + PAD + eps
    need(slow > 0, 'saddle above reference')
    dm = F(b['delta']) + PAD
    rm = F(b['rho']) + PAD
    cm = F(b['c']) - PAD
    need(cm > rm and cm * rm >= dm + rm * rm / 2, 'minimum correction/convex bridge')
    mm = cm - rm
    min_gap = dm * dm / (2 * mm)
    need(min_gap < F('.00001168099'), 'minimum reference energy gap')
    need(shigh + min_gap < F('.00570116'), 'global upper barrier display')
    # The deliberately tiny local box requires only ||H||<=1 for coupling.
    S, R = F(1, 100000), F(1, 1000)
    abox, mbox = astar - S - R, mstar - S - R
    beta = S + S * S / 2
    lateral = mbox * R * R / 2 - beta * R - S * S / 2
    need(abox > 0 and mbox > 0 and beta < mbox * R and lateral > 0, 'local saddle box')
    signed = set()
    points = []
    for v in joins['points']:
        label = v['label']
        need(label in ('plus', 'minus') and label not in signed, 'distinct point labels')
        signed.add(label)
        p = F(v['signed_projection_on_normalized_q'][0]) - PAD if label == 'plus' else -F(v['signed_projection_on_normalized_q'][1]) - PAD
        longitudinal = p - rho
        radius = F(v['distance_to_saddle_candidate'][1]) + PAD + rho
        endpoint_energy = F(v['relative_energy_to_reference'][1]) + PAD
        hermite_slack = slow - endpoint_energy - radius**3 / 6
        face_slack = R * longitudinal - S * radius
        need(longitudinal > S and face_slack > 0 and hermite_slack > 0, 'exact-saddle attachment: ' + label)
        points.append({'label': label, 'longitudinal_lower': str(longitudinal), 'distance_upper': str(radius), 'Hermite_margin': str(hermite_slack), 'face_margin': str(face_slack)})
    need(signed == {'plus', 'minus'}, 'both sides required')
    legs = []
    for name, leg, expected_nodes in [('return', ret, 8), ('annihilation', ann, 19)]:
        need(leg['nodes'] == expected_nodes and leg['segments'] == expected_nodes - 1 and leg['all_time_intervals_covered'], 'complete node/segment ledger')
        slack = slow - F(leg['path_relative_reference_upper']) - PAD
        need(slack > 0, 'leg crosses saddle level: ' + name)
        legs.append({'leg': name, 'strict_margin_below_saddle': str(slack)})
    need(ann['endpoint_charge_count'] == 0 and ann['monotone_continuation_to_constant'] and ann['compact_integer_lift_max'] == 0, 'terminal lift contract')
    paths = [pre / 'baseline.json', pre / 'saddle.json', con / 'join_points.json', con / 'leg_return.json', con / 'leg_annihilation.json']
    return {'status': 'CONDITIONAL_ANALYTIC_SCALAR_AUDIT_PASS',
            'scope': 'Exact scalar implications conditional on supplied local interval outputs; no independent residual or Hessian evaluation',
            'input_kind': 'archived or freshly supplied results, as named by result_root',
            'result_root': str(root), 'serialization_pad': str(PAD), 'Hessian_Lipschitz_majorant': '1',
            'source_hashes': {str(p.relative_to(root)): sha(p) for p in paths},
            'invariant_ball_slack': str(invariant_slack), 'saddle_inverse_gap_lower': str(min(astar, mstar)),
            'minimum_reference_gap_upper_rederived': str(min_gap), 'saddle_energy_error_upper_rederived': str(eps),
            'exact_saddle_relative_reference_lower': str(slow), 'exact_saddle_relative_reference_upper': str(shigh),
            'upper_barrier_rederived': str(shigh + min_gap), 'local_lateral_margin': str(lateral),
            'Hermite_attachments': points, 'global_legs': legs,
            'global_optimality_derived_here': False, 'heteroclinic_orbit_derived_here': False}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path, default=HERE / 'scalar_implications_archived.json')
    args = parser.parse_args()
    result = audit(args.results.resolve())
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_hashes', 'Hermite_attachments', 'global_legs')}, indent=2))


if __name__ == '__main__':
    main()
