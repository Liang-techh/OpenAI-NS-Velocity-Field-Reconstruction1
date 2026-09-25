"""Compare physical divergence and complete momentum after the slice lift."""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from radial_continuation import ROOT


def run(slice_filename='coupled_five_moment_slice.json',
        output_name='coupled_moment_physical_screen.json',
        sample_X=(1.005, 1.01, 1.1, 2., 2.99)):
    tau = .5*2**(-5.5)
    lift = CoupledMomentPhysicalLift(slice_filename=slice_filename)
    X, eta = np.meshgrid(sample_X,
                         [.2, .25, .3], indexing='ij')
    points = lift.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau)
    rows = []
    for label, field in (('before', lift.base), ('lifted', lift)):
        residual, divergence = independent_fd(
            field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        data = dict(field=label,
                    max_full_momentum=float(np.max(np.linalg.norm(
                        residual, axis=1))),
                    max_divergence=float(np.max(np.abs(divergence))),
                    nodes=[dict(X=float(x), eta=float(e),
                                residual_norm=float(np.linalg.norm(R)),
                                divergence=float(d))
                           for x, e, R, d in zip(X.ravel(), eta.ravel(),
                                                 residual, divergence)])
        rows.append(data)
        print(json.dumps({key: data[key] for key in
                          ('field', 'max_full_momentum',
                           'max_divergence')}), flush=True)
    report = dict(tau=tau, rows=rows,
                  slice_filename=slice_filename,
                  taper_width=lift.width, u_degree=lift.degree,
                  spatial_step=.001*np.sqrt(lift.nu*tau),
                  time_step=.00025*tau,
                  scope=f'{len(points)} physical finite-difference nodes in the fitted/interpolated axial band. No volume-L2, continuous cone, time-uniform or PDE acceptance. Pressure is inherited unchanged from the base field.',
                  accepted=False)
    (ROOT/output_name).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
