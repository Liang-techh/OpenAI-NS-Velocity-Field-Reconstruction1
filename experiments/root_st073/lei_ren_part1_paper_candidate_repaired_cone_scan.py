"""Same-source finite sampling of repaired Section 10.2 stress cones.

Ratios are formed in the retained pressure/width ring before evaluation.
This is a diagnostic, not an interval or whole-axis cone certificate.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_five_bump_inverse import build_problem

HERE = Path(__file__).resolve().parent


def nominal(value):
    return value.evaluate() if hasattr(value, 'evaluate') else value


def cone(field):
    if field['stress'] is None:
        raise ArithmeticError(field.get('stress_error'))
    F = field['F']
    stress = field['stress']
    # Retain small pressure/width atoms while dividing by the tiny F.
    st, sz, tt, tz = (stress[key] / F for key in
                      ('S_theta', 'S_z', 'T_theta', 'T_z'))
    kappa = -(st * st + sz * sz) / st
    dot = tt * st + tz * sz
    cross = -tt * sz + tz * st
    kn, dn = nominal(kappa), nominal(dot)
    prerequisites = nominal(F) > 0 and nominal(st) < 0
    if kn > 2:
        margin = 2 * dot * dot - (kappa - 2) * cross * cross
        branch = 'kappa>2'
    else:
        margin = -dot / (-st) - (2 - kappa)
        branch = 'kappa<=2'
    mn = nominal(margin)
    relaxed = prerequisites and dn < 0 and mn > 0
    def number(v):
        return mp.nstr(nominal(v), 65)
    return dict(kappa=number(kappa), kappa_minus_two=number(kappa - 2),
                normalized_S_theta=number(st), normalized_S_z=number(sz),
                normalized_T_theta=number(tt), normalized_T_z=number(tz),
                T_dot_S_over_F_squared=number(dot), branch=branch,
                branch_margin=number(margin), prerequisites_passed=bool(prerequisites),
                relaxed_passed=bool(relaxed),
                admissible_passed=bool(relaxed and kn > 2))


def run():
    source, centered, reference, moment_map, repaired = build_problem()
    rows = []
    with mp.workdps(moment_map.precision):
        points = {mp.mpf(s) for s in ('1', '1.1', '1.375', '1.625', '1.9', '2')}
        for center in moment_map.centers:
            for step in range(-5, 6):
                points.add(center + moment_map.radius * step / 5)
        for x in sorted(points):
            baseline = reference.evaluate_x(x, source.center)
            field = repaired.evaluate_x(x, source.center)
            if field['P0'] != baseline['P0'] or field['P0_Z'] != baseline['P0_Z']:
                raise AssertionError('Prescribed pressure datum changed')
            row = dict(x=mp.nstr(x, 40), baseline=cone(baseline), repaired=cone(field))
            rows.append(row)
            print('Repaired cone x', row['x'], 'kappa', row['repaired']['kappa'],
                  'relaxed', row['repaired']['relaxed_passed'], flush=True)
        failed = [row['x'] for row in rows if not row['repaired']['relaxed_passed']]
        inputs = ('lei_ren_part1_paper_candidate_five_bump_inverse.py',
                  'lei_ren_part1_paper_candidate_five_bump_inverse.json',
                  'lei_ren_part1_paper_candidate_transition_analytic.json',
                  'lei_ren_part1_paper_five_bump_inverse_field.py',
                  'lei_ren_part1_paper_five_bump_field.py',
                  'lei_ren_part1_paper_five_moment_reference_background.py',
                  'lei_ren_part1_paper_mp_stress.py',
                  'lei_ren_part1_paper_five_bump_map.py')
        report = dict(center='.3', Lambda='1e120', precision=moment_map.precision,
                      actual_candidate_used=True, samples=len(rows), rows=rows,
                      failed_relaxed_samples=failed,
                      all_sampled_relaxed_passed=not failed,
                      minimum_sampled_kappa=min(rows, key=lambda r: mp.mpf(r['repaired']['kappa']))['repaired']['kappa'],
                      normalized_ratios_formed_before_ring_evaluation=True,
                      finite_pressure_width_ring=True,
                      sampling_includes_support_edges_and_interiors=True,
                      P0_and_first_Z_preserved_at_samples=True,
                      uniform_cone_certified=False, source_error_enclosed=False,
                      quadrature_error_enclosed=False, whole_axis=False,
                      temporal_recursion=False, full_momentum_residual_validated=False,
                      input_hashes={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                    for name in inputs},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print('Same-candidate repaired cone scan complete:', len(rows), 'samples;',
              len(failed), 'relaxed failures; no uniform certificate', flush=True)
        return report


if __name__ == '__main__':
    run()
