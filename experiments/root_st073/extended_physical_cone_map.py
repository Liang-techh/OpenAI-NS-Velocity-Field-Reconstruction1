"""Map the actual full-residual stress cone for the extended shear bridge.

This differs from the paper-normalized leading-profile cone: the target is
the radial primitive of the finite-time *complete physical* residual.
Only a few nodes are sampled, so no continuous cone certificate follows.
"""
import json

import numpy as np

from high_frequency_shear_screen import make_field
from radial_peak_cone import operator, stress_primitive
from radial_continuation import ROOT


def physical_cone(field, X, eta, tau, order=16):
    point = field.compact.joined.inner.from_similarity(
        np.array([X]), np.array([eta]), tau)
    radius, _, z = point[0]
    u, J, R = operator(field, point, tau)
    target = stress_primitive(field, radius, z, tau, order=order)
    F = u[0, 1]/radius
    shear = np.array([J[0, 1, 0]-F, J[0, 2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    dot_n = float(target@N)
    dot_k = float(target@K)
    ratio = (float(abs(np.sqrt(lam2)/(2*F*N[0])*dot_k/dot_n))
             if lam2 > 0 and abs(2*F*N[0]*dot_n) > 1e-14 else None)
    return dict(X=X, eta=eta, order=order,
                lambda_squared=float(lam2),
                target_dot_N=dot_n, target_dot_K=dot_k,
                ratio=ratio,
                strict_pass=bool(lam2 > 0 and dot_n < 0 and
                                 ratio is not None and ratio < 1),
                residual_norm=float(np.linalg.norm(R[0])))


def run():
    tau = .5*2**(-5.5)
    field = make_field(16, 2.)
    rows = []
    for eta in (.2, .3):
        for X in (.75, .9, 1., 1.1, 1.25):
            row = physical_cone(field, X, eta, tau, order=64)
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, rows=rows,
                  scope='Ten-node Gauss64 physical full-residual cone map for extended N16 shear candidate. Not the normalized leading-profile theorem, a continuous cone, or a supported wave.',
                  accepted=False)
    (ROOT/'extended_physical_cone_map.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
