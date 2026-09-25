"""Screen a moment-null remote swirl shape against the physical stress analog."""
import json

from joined_field import ROOT
from moment_matched_joined_field import MomentMatchedJoinedField
from paper_moment_bridge import slice_data
from remote_moment_stress_screen import points_for_eta, summarize


def run():
    tau = .0084
    rows = []
    for amplitude in (-1., -.5, 0., .5, 1.):
        try:
            field = MomentMatchedJoinedField(null_amplitude=amplitude)
        except ValueError as exc:
            rows.append({'null_amplitude': amplitude, 'error': str(exc)})
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
        rows.append({'null_amplitude': amplitude, 'samples': samples,
                     'local_pass_count': sum(point['strict_local_pass']
                                             for row in samples
                                             for point in row['stress'])})
        print(json.dumps({'null_amplitude': amplitude,
                          'local_pass_count': rows[-1]['local_pass_count'],
                          'X10_eta0': samples[0]['stress'][1]['target_dot_N'],
                          'max_moment_defect': max(abs(row[key]) for row in samples
                                                   for key in ('angular_moment',
                                                               'kinetic_moment'))}),
              flush=True)
    report = {'rows': rows,
              'scope': 'Finite-amplitude moment-null remote swirl family, sampled physical stress-cone analog and direct slice moments. Not the normalized paper cone or a full residual/energy/force acceptance.',
              'accepted': False}
    (ROOT/'remote_null_shape_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
