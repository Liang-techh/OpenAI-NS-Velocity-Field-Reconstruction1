"""Check whether a locally passing wider-support candidate passes nearby."""
import json

import numpy as np

from compact_potential import CompactPotentialField
from curl_wave_cone_region import evaluate_height
from joined_field import ROOT
from radial_peak_cone import current_field


def run():
    field = current_field()
    tau = .5/64
    _, r_support, _, z_support = CompactPotentialField().support(tau)
    radial_axis = np.array([.0015, .002, .0025, .00275, .003,
                            .00325, .0035, .00375, .004, .0045])
    axial_axis = np.array([.0008, .0012, .0015, .002, .0025,
                           .003, .0035, .004])
    radial_profile = evaluate_height(field, radial_axis, .0025, tau,
                                     r_support, z_support)
    axial_profile = [evaluate_height(field, np.array([.0015]), z, tau,
                                     r_support, z_support)[0]
                     for z in axial_axis]
    candidate = json.loads((ROOT/'compact_potential'/'radial_peak_cone_r003_z0025.json').read_text())
    report = {'tau': tau,
              'candidate_point': candidate['point'],
              'candidate_local_cone': candidate['cone'],
              'candidate_covariance_nnls': candidate['covariance_nnls'],
              'radial_profile_z': .0025,
              'radial_profile': radial_profile,
              'axial_profile_r': .0015,
              'axial_profile': axial_profile,
              'scope': 'Physical local cone line screens only. A passing center with exact positive covariance does not imply a uniformly passing spatial envelope. Cone signs near tiny target values are sensitive to numerical error; no global support or NS acceptance follows.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_cone_patch_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'radial_passing': [x['r'] for x in radial_profile
                                           if x['strict_local_pass']],
                      'axial_passing': [x['z'] for x in axial_profile
                                          if x['strict_local_pass']],
                      'candidate_cone': candidate['cone'],
                      'candidate_covariance': candidate['covariance_nnls']}),
          flush=True)


if __name__ == '__main__':
    run()
