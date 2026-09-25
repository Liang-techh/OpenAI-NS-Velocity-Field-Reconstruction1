"""Physical radial cone span of narrow and plateau pressure candidates.

The twelve fit nodes stop at X=1.02. This screen integrates the complete
mean momentum primitive farther into the transition at eta=.25, retaining
both pressure gradients and all nonlinear ten-mode velocity terms.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_outer_radial_extension import make_extended_modes
from delayed_plateau_pressure_screen import (PlateauPressurePatch,
                                              RADIAL_SHAPES)
from delayed_similarity_curl_screen import CurlPatchedLift
from delayed_similarity_pressure_screen import (RADIAL_INTERVALS,
                                                SimilarityPressurePatch)
from delayed_swirl_cone_response import blocks, stress_primitive
from joint_collar_fit import kinematics
from radial_continuation import ROOT


XS = (1.008, 1.012, 1.016, 1.02, 1.03, 1.04, 1.05, 1.07,
      1.09, 1.12, 1.16, 1.2)
ETA = .25


def cone_row(u, grad, residual, block):
    center = block['center']
    target = stress_primitive(residual, block)
    F = u[center, 1]/block['radius']
    shear = np.array([grad[center, 1, 0]-F,
                      grad[center, 2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    dot_n = float(target@N)
    dot_k = float(target@K)
    if lam2 > 0 and abs(2*F*N[0]) > 1e-14:
        multiplier = np.sqrt(lam2)/(2*F*N[0])
        side = abs(multiplier*dot_k)
        strict_margin = -dot_n-side
        robust_margin = -.8*dot_n-side-.1
    else:
        strict_margin = robust_margin = -1e12
    return dict(X=block['X'], eta=block['eta'],
                radius=float(block['radius']),
                lambda_squared=float(lam2),
                target=target.tolist(),
                strict_margin=float(strict_margin),
                robust_margin=float(robust_margin),
                strict_pass=bool(strict_margin > 0),
                robust_pass=bool(robust_margin >= 0),
                residual_norm=float(np.linalg.norm(
                    residual[center])))


def run():
    mean_source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    pressure_source = json.loads((ROOT/'delayed_outer_pressure_fit.json').read_text())
    plateau_source = json.loads((ROOT/'delayed_plateau_pressure_screen.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = CurlPatchedLift(base, make_extended_modes(base),
                           mean_source['coefficients'])
    old = SimilarityPressurePatch(base, np.zeros(6))
    plateau = PlateauPressurePatch(base)
    tau = .5*2**(-5.5)
    radial_edges = [edge for interval in RADIAL_INTERVALS
                    for edge in interval]
    radial_edges += [edge for rise, fall in RADIAL_SHAPES
                     for edge in (*rise, *fall)]
    radial_edges += [1.01, 1.05]
    points, segments = blocks(base, tau, order=32, xs=XS,
                              etas=(ETA,),
                              extra_radial_edges=radial_edges)
    hs, ht = .0005*np.sqrt(base.nu*tau), .0001*tau
    u, grad, part = kinematics(mean, points, tau, hs, ht)
    mean_residual = part+np.einsum('pab,pb->pa', grad, u)
    _, old_gradients = old.basis(points, tau)
    _, plateau_gradients = plateau.basis(points, tau)
    pressure_gradients = np.concatenate((old_gradients,
                                         plateau_gradients), axis=2)
    old_coefficients = np.r_[pressure_source['pressure_coefficients'],
                              np.zeros(plateau.count)]
    cap_coefficients = np.asarray(
        plateau_source['searches']['combined_cap_1e3']['coefficients'])
    candidates = dict(mean_only=np.zeros(6+plateau.count),
                      old_six=old_coefficients,
                      plateau_cap_1e3=cap_coefficients)
    reports = {}
    for name, coefficients in candidates.items():
        residual = mean_residual+pressure_gradients@coefficients
        rows = [cone_row(u, grad, residual, block)
                for block in segments]
        reports[name] = rows
        print(json.dumps(dict(name=name,
                              robust_pass=[row['X'] for row in rows
                                           if row['robust_pass']],
                              strict_pass=[row['X'] for row in rows
                                           if row['strict_pass']],
                              max_center_residual=max(
                                  row['residual_norm'] for row in rows))),
              flush=True)
    report = dict(source_mean='delayed_outer_admission_search.json',
                  source_pressure='delayed_outer_pressure_fit.json',
                  source_plateau='delayed_plateau_pressure_screen.json',
                  tau=tau, eta=ETA, X=XS, quadrature_order=32,
                  candidates=reports,
                  scope='One-eta radial screen with direct full '
                        'Cartesian momentum primitive and analytic '
                        'pressure gradients. Extra radii were not used '
                        'in the pressure LP. Coarse finite nodes only: '
                        'no continuous cone, supported wave, moment '
                        'closure, or PDE acceptance.', accepted=False)
    (ROOT/'delayed_plateau_cone_span.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
