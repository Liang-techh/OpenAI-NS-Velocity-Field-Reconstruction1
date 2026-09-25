"""Materialize the screened fixed-moment momentum-tangent profile."""

import argparse
import copy
import json

from radial_continuation import ROOT


def run(screen_name='delayed_momentum_tangent_screen.json',
        output_name='delayed005_wide04_momentum_tangent.json'):
    source_name = 'delayed005_wide04_reoptimized_E.json'
    source = json.loads((ROOT/source_name).read_text())
    screen = json.loads((ROOT/screen_name).read_text())
    if screen['slice_filename'] != source_name:
        raise ValueError('Momentum screen used a different slice source')
    candidate = copy.deepcopy(source)
    row = next(row for row in candidate['rows'] if row['eta'] == .3)
    coefficients = screen['selected']['coefficients']
    if len(coefficients) != len(row['u_coefficients']):
        raise ValueError('Momentum tangent changed profile degree')
    row['u_coefficients'] = coefficients
    row['variant'] = 'momentum_tangent'
    for stale in ('optimized_curvature_proxy', 'curvature_ratio',
                  'optimizer_success', 'final_delta'):
        row.pop(stale, None)
    original_row = next(row for row in source['rows'] if row['eta'] == .3)
    defect = screen.get('selected_slice_moment_defect',
                        screen['selected'].get('slice_moment_defect'))
    if defect is None:
        raise ValueError('Momentum screen lacks a slice-moment check')
    row['max_abs_five_moment_defect'] = (
        original_row['max_abs_five_moment_defect']
        + max(abs(value) for value in defect))
    row['moment_defect_method'] = 'coefficient quadrature; physical check is separate'
    row['five_moments_restored'] = True
    row['momentum_tangent_screen'] = screen_name
    candidate['scope'] = (
        'Two fixed-eta five-moment slices at the reference time; eta=.3 '
        'uses a projected momentum-tangent U correction. Local Cartesian '
        'screens are separate and no global cone or PDE acceptance is claimed.')
    candidate['source_slice'] = source_name
    candidate['momentum_tangent_screen'] = screen_name
    candidate['accepted'] = False
    (ROOT/output_name).write_bytes(
        (json.dumps(candidate, indent=2)+'\n').encode())
    print(output_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--screen-name',
                        default='delayed_momentum_tangent_screen.json')
    parser.add_argument('--output-name',
                        default='delayed005_wide04_momentum_tangent.json')
    args = parser.parse_args()
    run(args.screen_name, args.output_name)
