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
                a = -st/f if f else float('nan')
                t0 = -sz/st if st else float('nan')
                p1 = float(value.get('I_theta', 0.))/f if f else float('nan')
                p2 = float(value.get('I_z', 0.))/f if f else float('nan')
                H = p1+p2*t0
                tt, tz = float(value.get('T_theta',0.)), float(value.get('T_z',0.))
                dot = tt*st+tz*sz
                perpendicular = -tt*sz+tz*st
                kappa = value['kappa']
                relaxed_margin = (2*dot*dot-(float(kappa)-2)*perpendicular**2
                                  if kappa is not None and kappa>2
                                  else -dot-(-f*st)*(2-float(kappa or 0.)))
                rows.append({'R': r, 'Z': value['Z'], 'F': f,
                             'S_theta': st, 'S_z': sz, 'kappa': value['kappa'],
                             'minimum_abs_U_R_boundary': (max(0., required_sq)/(2*r))**.5,
                             'actual_abs_U_R': abs(sz)/(2*r)**.5,
                             'a': a, 't0': t0, 'H_t0': H,
                             'frozen_loop_necessary_margin_H_minus_2': H-2,
                             'frozen_loop_necessary_input_pass': bool(a>0 and H>2),
                             'actual_relaxed_cone_margin': relaxed_margin,
                             'actual_relaxed_cone_pass': bool(f>0 and st<0 and dot<0 and relaxed_margin>0),
                             'actual_stress_component_max': max(abs(tt),abs(tz)),
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
              'original_candidate_metadata_unrecorded': data.get('original_candidate_metadata_unrecorded', False),
              'input_candidate_metadata': data.get('parameters', {}).get('matched_metadata'),
              'necessary_condition': 'S_z^2 > -2 F S_theta - S_theta^2; also require F>0 and S_theta<0',
              'weak_anchor_obstruction': 'On F=G(Z) R^(-p), U_R=0 implies kappa=2p. For p=.005 this is .01, below 2, regardless of terminal moment closure.',
              'next_action': 'Restore axial radial shear across weak-anchor gaps while preserving mass, mixed, quadratic and angular moments; then recompute actual stress direction.',
              'source_shear_loop_gate': 'Lei-Ren Section 11 additionally requires relaxed-cone input with H(t0)>2 and positive boundary margins. Raising kappa alone is insufficient.',
              'rows': rows,
              'necessary_shear_failed_count': sum(row['necessary_shear_failed'] for row in rows),
              'frozen_loop_necessary_input_failed_count': sum(not row['frozen_loop_necessary_input_pass'] for row in rows),
              'actual_relaxed_cone_failed_count': sum(not row['actual_relaxed_cone_pass'] for row in rows),
              'scope': 'Diagnostic on measured points with nonzero stress; zero-stress core requires a separate directional boundary treatment. Neither scalar screen establishes the relaxed or admissible cone.',
              'whole_cone_validated': False, 'pde_validated': False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','input_candidate_metadata')}))
    return report


if __name__ == '__main__': run()
