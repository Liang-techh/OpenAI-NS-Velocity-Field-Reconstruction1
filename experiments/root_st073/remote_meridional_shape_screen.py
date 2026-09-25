"""Screen a second remote streamfunction shape under both radial moments."""
import json

from joined_field import ROOT
from moment_matched_joined_field import MomentMatchedJoinedField
from paper_moment_bridge import slice_data
from remote_moment_stress_screen import points_for_eta, summarize


def run():
    tau = .0084
    rows = []
    for amplitude in (-5., -3., -1., 0., 1., 3.):
        try:
            field = MomentMatchedJoinedField(
                meridional_null_amplitude=amplitude)
        except ValueError as exc:
            rows.append({'meridional_null_amplitude': amplitude,
                         'error': str(exc)})
            continue
        samples = []
        for eta in (0., .2, .35):
            parts, points = points_for_eta(field, eta, tau)
            stress, maximum = summarize(field, parts, points, tau)
            moment = slice_data(field, eta, tau, 48,
                                field.patch_start, field.patch_end)
            samples.append({'eta': eta, 'stress': stress,
                            'sample_max_residual': maximum,
                            'angular_moment': moment['baseline_angular_moment'],
                            'kinetic_moment': moment['kinetic_moment_baseline']})
        row = {'meridional_null_amplitude': amplitude,
               'samples': samples,
               'local_pass_count': sum(point['strict_local_pass']
                                       for sample in samples
                                       for point in sample['stress'])}
        rows.append(row)
        print(json.dumps({'meridional_null_amplitude': amplitude,
                          'local_pass_count': row['local_pass_count'],
                          'X10_eta0': samples[0]['stress'][1]['target_dot_N'],
                          'max_moment_defect': max(abs(sample[key])
                                                   for sample in samples
                                                   for key in ('angular_moment',
                                                               'kinetic_moment'))}),
              flush=True)
    report = {'rows': rows,
              'scope': 'Six-point meridional streamfunction family with physical stress-cone analog at 9 points per field. Direct slice moments but no continuum, normalized paper cone, wave, force or energy acceptance.',
              'accepted': False}
    (ROOT/'remote_meridional_shape_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
