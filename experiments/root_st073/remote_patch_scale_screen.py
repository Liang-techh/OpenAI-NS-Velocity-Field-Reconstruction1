"""Move a two-moment patch outward and evaluate the normalized relaxed cone."""
import json

from joined_field import ROOT
from moment_matched_joined_field import MomentMatchedJoinedField
from normalized_relaxed_cone import cone_point
from paper_moment_bridge import slice_data


def run():
    tau = .0084
    rows = []
    for patch in (None, (16., 64.), (64., 256.),
                  (256., 1024.), (1024., 4096.)):
        field = MomentMatchedJoinedField(patch=patch)
        start, end = field.patch_start, field.patch_end
        radii = [start+(end-start)*fraction for fraction in
                 (1/6, .5, 5/6)]
        points = []
        for eta in (0., .2, .35):
            moment = slice_data(field, eta, tau, 48, start, end)
            for X in radii:
                point = cone_point(field, X, eta, tau)
                point['angular_moment'] = moment['baseline_angular_moment']
                point['kinetic_moment'] = moment['kinetic_moment_baseline']
                points.append(point)
        row = {'patch_X': [start, end], 'points': points,
               'relaxed_pass_count': sum(p['relaxed_pass'] for p in points),
               'admissible_pass_count': sum(p['admissible_pass'] for p in points),
               'sample_max_moment_defect': max(
                   abs(p[key]) for p in points for key in
                   ('angular_moment', 'kinetic_moment'))}
        rows.append(row)
        print(json.dumps({'patch_X': row['patch_X'],
                          'relaxed_pass_count': row['relaxed_pass_count'],
                          'admissible_pass_count': row['admissible_pass_count'],
                          'X_mid_eta0_vs': points[1].get('vs'),
                          'X_mid_eta0_Pc': points[1].get('Pc'),
                          'sample_max_moment_defect':
                              row['sample_max_moment_defect']}), flush=True)
    report = {'tau': tau, 'rows': rows,
              'scope': 'Five remote patch scales on finite-window normalized snapshots, with two radial moment diagnostics and paper relaxed-cone samples. The pressure is reconstructed only for the leading cone diagnostic, not installed into the physical field. No five-moment interface, continuum, finite-energy, force or full momentum acceptance.',
              'accepted': False}
    (ROOT/'remote_patch_scale_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
