"""Test one shifted solenoidal swirl mode on an expanded cone window.

The amplitude is selected from the shear-sign screen, then the full
physical momentum primitive is recomputed on five radii and three axial
slices. Existing plus plateau pressure modes are fitted by a constrained
component-minimax LP. This remains a sampled necessary-stage experiment.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import SimilaritySwirlMode, mode_profile
from delayed_momentum_tangent_screen import nodes, stats
from delayed_multimode_cone_fit import BASE_NAME, load_cache
from delayed_outer_radial_extension import (CACHE_NAME,
                                             make_extended_modes, signature)
from delayed_plateau_pressure_screen import (PlateauPressurePatch,
                                              RADIAL_SHAPES,
                                              minmax_pressure,
                                              plateau_stress_columns)
from delayed_pressure_cone_window import pressure_stress_columns
from delayed_similarity_curl_screen import CurlPatchedLift
from delayed_similarity_pressure_screen import (RADIAL_INTERVALS,
                                                SimilarityPressurePatch)
from delayed_swirl_cone_response import blocks, stress_primitive
from joint_collar_fit import kinematics
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


XS = (1.008, 1.012, 1.016, 1.02, 1.03)
ETAS = (.2, .25, .3)
NEW_RADIAL = (1.02, 1.06)
NEW_AXIAL = (.14, .36)
AMPLITUDE = -.0105


def geometry(u, grad, radius):
    F = u[1]/radius
    shear = np.array([grad[1, 0]-F, grad[2, 0]])
    N = shear/np.linalg.norm(shear)
    K = np.array([-N[1], N[0]])
    lam2 = -2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
    multiplier = (np.sqrt(lam2)/(2*F*N[0])
                  if lam2 > 0 and abs(2*F*N[0]) > 1e-14 else None)
    return N, K, float(lam2), multiplier


def run():
    mean_source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    original_modes = make_extended_modes(base)
    extra_mode = SimilaritySwirlMode(base, NEW_RADIAL,
                                     eta_interval=NEW_AXIAL)
    tau = .5*2**(-5.5)
    mean = CurlPatchedLift(base, original_modes+[extra_mode],
                           [*mean_source['coefficients'], AMPLITUDE])
    old_pressure = SimilarityPressurePatch(base, np.zeros(6))
    plateau = PlateauPressurePatch(base)
    radial_edges = [edge for interval in RADIAL_INTERVALS
                    for edge in interval]
    radial_edges += [edge for rise, fall in RADIAL_SHAPES
                     for edge in (*rise, *fall)]
    radial_edges += [1.01, 1.05, *NEW_RADIAL]
    points, segments = blocks(base, tau, order=32, xs=XS,
                              etas=ETAS,
                              extra_radial_edges=radial_edges)
    hs, ht = .0005*np.sqrt(base.nu*tau), .0001*tau
    u, grad, part = kinematics(mean, points, tau, hs, ht)
    residual = part+np.einsum('pab,pb->pa', grad, u)
    centers = np.array([block['center'] for block in segments])
    center_points = points[centers]
    center_residual = residual[centers]
    A, rhs, geometry_rows = [], [], []
    for block in segments:
        i = block['center']
        target = stress_primitive(residual, block)
        N, K, lam2, multiplier = geometry(
            u[i], grad[i], block['radius'])
        if lam2 <= 0 or multiplier is None:
            raise ValueError(f'Invalid lambda at {block["X"]},{block["eta"]}')
        old_columns = pressure_stress_columns(
            old_pressure, points[i], tau)
        new_columns = plateau_stress_columns(
            plateau, points[i], tau)
        columns = np.column_stack((old_columns, new_columns))
        for sign in (-1., 1.):
            direction = .8*N+sign*multiplier*K
            A.append(direction@columns)
            rhs.append(-.1-direction@target)
        geometry_rows.append(dict(X=block['X'], eta=block['eta'],
                                  lambda_squared=lam2,
                                  target=target.tolist(),
                                  residual_norm=float(np.linalg.norm(
                                      residual[i]))))
    A, rhs = np.asarray(A), np.asarray(rhs)
    nearby_tau = .5*2**(-5.4)
    nearby_points, _, _ = nodes(base, (1.01, 1.015, 1.025),
                                 (.22, .28, .32), nearby_tau)
    nearby_hs, nearby_ht = (.001*np.sqrt(base.nu*nearby_tau),
                            .00025*nearby_tau)
    nearby_u, nearby_g, nearby_part = kinematics(
        mean, nearby_points, nearby_tau, nearby_hs, nearby_ht)
    nearby_residual = (nearby_part+
                       np.einsum('pab,pb->pa', nearby_g, nearby_u))
    _, old_train_g = old_pressure.basis(center_points, tau)
    _, new_train_g = plateau.basis(center_points, tau)
    _, old_near_g = old_pressure.basis(nearby_points, nearby_tau)
    _, new_near_g = plateau.basis(nearby_points, nearby_tau)
    train_gradient = np.concatenate((old_train_g, new_train_g), axis=2)
    near_gradient = np.concatenate((old_near_g, new_near_g), axis=2)
    gradient_matrix = np.concatenate((
        train_gradient.reshape(-1, 18),
        near_gradient.reshape(-1, 18)))
    vector = np.r_[center_residual.ravel(), nearby_residual.ravel()]
    pressure_fits = {}
    for name, active, limit in (
            ('old_six', list(range(6)), None),
            ('all_eighteen', list(range(18)), None),
            ('all_cap_1e3', list(range(18)), 1e3)):
        fit = minmax_pressure(A, rhs, gradient_matrix,
                              vector, active, limit)
        if fit['feasible']:
            pc = np.asarray(fit['coefficients'])
            fit.update(train_max=stats(center_residual+
                                       train_gradient@pc)['max'],
                       nearby_max=stats(nearby_residual+
                                        near_gradient@pc)['max'])
        pressure_fits[name] = fit
        print(json.dumps(dict(name=name,
                              feasible=fit['feasible'],
                              component_max=fit.get('actual_component_max'),
                              train_max=fit.get('train_max'),
                              nearby_max=fit.get('nearby_max'))), flush=True)
    _, moments = load_cache(ROOT/CACHE_NAME, signature())
    coeff = np.asarray(mean_source['coefficients'])
    moment_data = []
    for row in moments:
        du, de = mode_profile(extra_mode, base, row['X'],
                              row['eta'], tau)
        U = row['U']+row['Bu']@coeff+AMPLITUDE*du
        E = row['E']+row['Be']@coeff+AMPLITUDE*de
        defect = moment_vector(U, E, row['X'], row['weights'])-row['target']
        moment_data.append(dict(eta=row['eta'],
                                defect=defect.tolist(),
                                min_relative_E=float(np.min(E/row['E0']))))
    report = dict(source='delayed_outer_admission_search.json',
                  new_radial=NEW_RADIAL, new_axial=NEW_AXIAL,
                  amplitude=AMPLITUDE, tau=tau, X=XS, eta=ETAS,
                  quadrature_order=32,
                  geometry=geometry_rows,
                  train_before=stats(center_residual),
                  nearby_before=stats(nearby_residual),
                  moments=moment_data,
                  pressure_fits=pressure_fits,
                  scope='One shifted exact-solenoidal swirl mode with '
                        'full physical stress primitive on an expanded '
                        '15-node cone, outgoing moment audit, and '
                        'sampled pressure-constrained momentum LP. '
                        'No continuous cone, wave, volume L2, or PDE '
                        'acceptance.', accepted=False)
    (ROOT/'delayed_shifted_swirl_cone_fit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(lambda_X103=[
        row['lambda_squared'] for row in geometry_rows
        if row['X'] == 1.03],
        train_before=report['train_before'],
        nearby_before=report['nearby_before'],
        largest_moment_defect=max(abs(value)
                                  for row in moment_data
                                  for value in row['defect']))), flush=True)


if __name__ == '__main__':
    run()
