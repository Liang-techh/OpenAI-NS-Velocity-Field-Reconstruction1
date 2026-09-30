"""Necessary shear screen of the measured candidate; no cone certificate."""
import json
from pathlib import Path


def run():
    source = Path(__file__).with_name('lei_ren_part1_extended_profile_stress_checks.json')
    data = json.loads(source.read_text(encoding='utf-8'))
    # With S_theta=2R F_R and S_z=sqrt(2R) U_R,
    # kappa>2 requires S_z^2 > -2F S_theta-S_theta^2.
    # This is necessary only; stress direction imposes further conditions.
    rows = []
    def visit(value):
        if isinstance(value, dict):
            if all(key in value for key in ('R','Z','F','S_theta','S_z','kappa')):
                r, f, st, sz = (float(value[k]) for k in ('R','F','S_theta','S_z'))
                required_sq = -2*f*st-st*st
                rows.append({'R': r, 'Z': value['Z'], 'F': f,
                             'S_theta': st, 'S_z': sz, 'kappa': value['kappa'],
                             'minimum_abs_U_R_boundary': (max(0., required_sq)/(2*r))**.5,
                             'actual_abs_U_R': abs(sz)/(2*r)**.5,
                             'necessary_shear_failed': bool(required_sq>=sz*sz)})
            else:
                for child in value.values(): visit(child)
        elif isinstance(value, list):
            for child in value: visit(child)
    # Use only the main grid, rather than duplicating the heat-tail rows.
    for key, value in data.items():
        if key not in ('tail_diagnostic', 'heat_tail_diagnostic'): visit(value)
    unique = {(row['R'], row['Z']): row for row in rows}
    rows = list(unique.values())
    if not rows: raise ValueError('No measured stress rows found')
    report = {'input_receipt': source.name,
              'original_candidate_metadata_unrecorded': True,
              'necessary_condition': 'S_z^2 > -2 F S_theta - S_theta^2; also require F>0 and S_theta<0',
              'weak_anchor_obstruction': 'On F=G(Z) R^(-p), U_R=0 implies kappa=2p. For p=.005 this is .01, below 2, regardless of terminal moment closure.',
              'next_action': 'Restore axial radial shear across weak-anchor gaps while preserving mass, mixed, quadratic and angular moments; then recompute actual stress direction.',
              'rows': rows,
              'necessary_shear_failed_count': sum(row['necessary_shear_failed'] for row in rows),
              'whole_cone_validated': False, 'pde_validated': False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
    return report


if __name__ == '__main__': run()
