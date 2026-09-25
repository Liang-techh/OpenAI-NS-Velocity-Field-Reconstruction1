"""Compare smooth U basis capacity as its onset enters the cone region.

Earlier onsets are diagnostic alternatives, not admissible additions to
the existing fifteen-node stress cone until that cone is re-audited.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
from delayed_remote_positive_e import quadrature
from delayed_remote_u_basis_capacity import capacity, u_basis
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


def run():
    source = json.loads((ROOT/'delayed_remote_swirl_floor_capacity.json').read_text())
    e_row = next(entry['result'] for entry in source['rows']
                 if entry['relative_E_floor'] == .05)
    e_source = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target_field = make_field(16, 2.)
    X, weights = quadrature()
    eta = .3
    tau = .5*2**(-5.5)
    U, E = profile(mean, X, eta, tau)
    U0, E0 = profile(target_field, X, eta, tau)
    target = moment_vector(U0, E0, X, weights)
    radial = np.array([bump(X, *interval)[0]
                       for interval in e_source['radial_intervals']])
    axial = bump(np.array([eta]), *e_source['mode_eta_interval'])[0][0]
    E_new = E+axial*np.asarray(e_row['coefficients'])@radial
    rows = []
    for start in (1.0, 1.01, 1.02, 1.025, 1.03):
        for width in (.005, .01, .02, .04):
            for degree in (11, 19, 31):
                B = u_basis(X, start, 3.5, width, degree)
                rows.append(dict(start=start, width=width,
                                 degree=degree,
                                 protected_cone_unchanged=(start>=1.03),
                                 **capacity(U, E_new, target, X,
                                            weights, B)))
    report = dict(source='delayed_remote_swirl_floor_capacity.json',
                  eta=eta, relative_E_floor=.05, end=3.5, rows=rows,
                  scope='Fixed-slice M/J KKT and S capacity for tapered '
                        'polynomial U corrections at varied onset. '
                        'Any start below 1.03 perturbs the existing '
                        'sampled inner stress cone and requires full '
                        'cone and momentum re-evaluation.',
                  accepted=False)
    (ROOT/'delayed_remote_u_start_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps([(s, min(r['degree'] for r in rows
                              if r['start']==s and r['width']==w
                              and r['target_minus_S_min']>=0)
                       if any(r['start']==s and r['width']==w
                              and r['target_minus_S_min']>=0 for r in rows)
                       else None, w)
                      for s in (1.0,1.01,1.02,1.025,1.03)
                      for w in (.005,.01,.02,.04)]), flush=True)


if __name__ == '__main__':
    run()
