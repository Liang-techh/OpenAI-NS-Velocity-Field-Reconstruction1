"""Screen two independent meridional streamfunction modes in the remote patch."""
import json

from joined_field import ROOT
from moment_matched_joined_field import MomentMatchedJoinedField
from paper_moment_bridge import slice_data
from remote_moment_stress_screen import points_for_eta, summarize


PAIRS = ((-5., -12.), (-5., 12.), (-3., -12.), (-3., 12.))


def run():
    tau = .0084
    rows = []
    for odd, even in PAIRS:
        try:
            field = MomentMatchedJoinedField(
                meridional_null_amplitude=odd,
                meridional_even_amplitude=even)
        except ValueError as exc:
            rows.append({'odd_amplitude': odd, 'even_amplitude': even,
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
        row = {'odd_amplitude': odd, 'even_amplitude': even,
               'samples': samples,
               'local_pass_count': sum(point['strict_local_pass']
                                       for sample in samples
                                       for point in sample['stress'])}
        rows.append(row)
        print(json.dumps({'odd_amplitude': odd, 'even_amplitude': even,
                          'local_pass_count': row['local_pass_count'],
                          'X10_eta0': samples[0]['stress'][1]['target_dot_N'],
                          'X10_eta0_ratio': samples[0]['stress'][1]['cone_ratio']}),
              flush=True)
    report = {'rows': rows,
              'scope': 'Four two-shape meridional fields, each at 9 physical stress-cone analog points. No continuum or normalized paper cone, wave, force or energy acceptance.',
              'accepted': False}
    (ROOT/'remote_two_shape_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
