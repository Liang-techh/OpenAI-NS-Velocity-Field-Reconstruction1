"""Check E positivity across the physical U-lift axial support."""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
from delayed_remote_positive_e import quadrature
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from paper_moment_bridge import bump
from radial_continuation import ROOT


def run():
    source = json.loads((ROOT/'delayed_remote_u_physical.json').read_text())
    e_source = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16,2.)
    X, _ = quadrature(order=48)
    radial = np.array([bump(X,*interval)[0]
                       for interval in e_source['radial_intervals']])
    c = np.asarray(source['slice']['e_coefficients'])
    tau = .5*2**(-5.5)
    rows = []
    for eta in (.26,.27,.28,.29,.3,.31,.32,.33,.34):
        _, E = profile(mean,X,eta,tau)
        _, E0 = profile(target,X,eta,tau)
        f = bump(np.array([eta]), *source['axial_interval'])[0][0]
        corrected = E+f*c@radial
        rows.append(dict(eta=eta,min_E=float(np.min(corrected)),
                         min_relative_E=float(np.min(corrected/E0))))
    report = dict(source='delayed_remote_u_physical.json',
                  eta_rows=rows,
                  scope='Nine eta slices on an independent Gauss-48 '
                        'grid. Positive samples do not certify '
                        'continuum positivity or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_remote_u_positivity_audit.json').write_bytes(
        (json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(rows),flush=True)


if __name__ == '__main__':
    run()
