"""Widen the fixed-start remote moment patch and screen the relaxed cone."""
import json

from joined_field import ROOT
from moment_matched_joined_field import MomentMatchedJoinedField
from normalized_relaxed_cone import cone_point
from paper_moment_bridge import slice_data


def run():
    tau = .0084
    rows = []
    for end in (64., 256., 4096., 16384.):
        field = MomentMatchedJoinedField(patch=(4., end))
        radii = [4.+(end-4.)*fraction for fraction in (1/6, .5, 5/6)]
        points = []
        for eta in (0., .2, .35):
            moment = slice_data(field, eta, tau, 48, 4., end)
            for X in radii:
                point = cone_point(field, X, eta, tau)
                point['angular_moment'] = moment['baseline_angular_moment']
                point['kinetic_moment'] = moment['kinetic_moment_baseline']
                points.append(point)
        row = {'patch_X': [4., end], 'points': points,
               'relaxed_pass_count': sum(p['relaxed_pass'] for p in points),
               'admissible_pass_count': sum(p['admissible_pass'] for p in points),
               'sample_max_moment_defect': max(abs(p[key]) for p in points
                                               for key in ('angular_moment',
                                                           'kinetic_moment'))}
        rows.append(row)
        print(json.dumps({'patch_X': row['patch_X'],
                          'relaxed_pass_count': row['relaxed_pass_count'],
                          'mid_eta0_vs': points[1].get('vs'),
                          'mid_eta0_upper': points[1].get('upper'),
                          'sample_max_moment_defect':
                              row['sample_max_moment_defect']}), flush=True)
    report = {'tau': tau, 'rows': rows,
              'scope': 'Four widening fixed-start remote patches; normalized relaxed-cone samples and two radial moment defects. Pi is reconstructed only for the leading diagnostic. No five-moment interface, full momentum, finite energy, force or continuum acceptance.',
              'accepted': False}
    (ROOT/'remote_patch_width_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
