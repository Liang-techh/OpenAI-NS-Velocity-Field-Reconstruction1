"""Locate scale-compatible momentum error before expanding spatial controls."""
import hashlib
import json
from pathlib import Path
import numpy as np
from wave_residual_harmonics import budget

ROOT = Path(__file__).resolve().parent


def run():
    cache = ROOT / 'scale_generator_momentum_defect.npz'
    source = ROOT / 'localized_drift_1800_fit.json'
    candidate = json.loads(source.read_text())
    defect_path = ROOT / 'scale_generator_momentum_defect.json'
    defect = json.loads(defect_path.read_text())
    if defect['status'] != 'completed':
        raise ValueError('Defect calculation must be complete')
    nu = defect['inputs']['viscosity']
    with np.load(cache) as data:
        points, weights = data['points'], data['weights']
        residual = data['generator_residual']
        velocity, jacobian = data['velocity'], data['gradient']
        terms = dict(generator_ut=data['generator_ut'],
                     convection=np.einsum('nij,nj->ni', data['gradient'], data['velocity']),
                     viscosity=-nu*data['laplacian'],
                     pressure=data['pressure_gradient'])
    radii = np.hypot(points[:, 0], points[:, 1])
    masks = []
    for name in ('first', 'second'):
        c = candidate['inputs'][name + '_patch_center']
        w = candidate['inputs'][name + '_patch_widths']
        masks.append((abs(radii-c[0]) < w[0]) & (abs(points[:, 2]-c[1]) < w[1]))
    regions = {'first_only': masks[0] & ~masks[1],
               'second_only': masks[1] & ~masks[0],
               'overlap': masks[0] & masks[1],
               'outside_both': ~(masks[0] | masks[1])}
    squares = weights * np.sum(residual**2, axis=1)
    total = float(squares.sum())
    rows = {}
    for name, mask in regions.items():
        rows[name] = dict(point_count=int(mask.sum()),
                          squared_L2_fraction=float(squares[mask].sum()/total),
                          volume_L2=float(np.sqrt(squares[mask].sum())))
    report = dict(status='completed', accepted=False, scale_recursion_established=False,
                  cache_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
                  candidate_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  defect_report_sha256=hashlib.sha256(defect_path.read_bytes()).hexdigest(),
                  angular=budget(points, weights, residual, 12), regions=rows,
                  scope='Frozen inner-grid error budget. Regions are disjoint supports, not claims of basis approximation capacity or PDE accuracy.')
    mean_swirl_square = report['angular']['modes'][0]['cylindrical_component_squared_L2'][1]
    report['axisymmetric_azimuthal_defect'] = dict(
        volume_L2=float(np.sqrt(mean_swirl_square)),
        squared_L2_fraction=float(mean_swirl_square/total),
        interpretation='The true angular mean of (1/r)*dtheta(p) vanishes for single-valued periodic pressure. This resolved-grid m=0 azimuthal residual therefore identifies a pressure-invariant momentum component. Numerical value is a sampled estimate, not a continuum lower-bound certificate.')
    theta = np.arctan2(points[:, 1], points[:, 0])
    def mean_theta(values):
        return (-np.sin(theta)*values[:, 0]+np.cos(theta)*values[:, 1]).reshape(-1,12).mean(axis=1)
    target = mean_theta(residual)
    ring_weights = weights.reshape(-1,12).sum(axis=1)
    report['axisymmetric_azimuthal_terms'] = {
        name: dict(volume_L2=float(np.sqrt(np.sum(ring_weights*mean_theta(values)**2))),
                   signed_projection_fraction=float(np.sum(ring_weights*mean_theta(values)*target)/mean_swirl_square))
        for name, values in terms.items()}
    report['axisymmetric_azimuthal_terms']['scope'] = 'Signed projections onto the mean azimuthal residual sum to one; term norms do not add. Viscosity read from the frozen defect source.'
    er = np.column_stack((np.cos(theta), np.sin(theta), np.zeros(len(theta))))
    et = np.column_stack((-np.sin(theta), np.cos(theta), np.zeros(len(theta))))
    def ring_mean(values):
        return values.reshape(-1,12).mean(axis=1)
    ur = ring_mean(np.sum(er*velocity, axis=1))
    uth = mean_theta(velocity)
    uz = ring_mean(velocity[:, 2])
    duth_dr = ring_mean(np.einsum('ni,nij,nj->n', et, jacobian, er))
    duth_dz = ring_mean(np.einsum('ni,ni->n', et, jacobian[:, :, 2]))
    mean_flow = ur*duth_dr + uz*duth_dz + ur*uth/radii.reshape(-1,12)[:,0]
    fluctuations = mean_theta(terms['convection']) - mean_flow
    report['mean_azimuthal_convection_split'] = {
        name: dict(volume_L2=float(np.sqrt(np.sum(ring_weights*values**2))),
                   signed_projection_onto_residual=float(np.sum(ring_weights*values*target)/mean_swirl_square))
        for name, values in [('convection_of_angular_mean',mean_flow),
                             ('fluctuation_momentum_flux',fluctuations)]}
    report['mean_azimuthal_convection_split']['scope'] = 'Angular Reynolds decomposition using resolved cylindrical means and Cartesian Jacobians. Fluctuation term is mean total convection minus convection of the angular mean. No sign/positivity or continuum certificate follows from term norms.'
    (ROOT / 'scale_generator_defect_budget.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(regions=rows, modes=report['angular']['modes'])))


if __name__ == '__main__':
    run()
