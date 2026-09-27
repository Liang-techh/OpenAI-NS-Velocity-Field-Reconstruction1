"""Fixed physical-cylinder observables for saved local state trajectories.

Sampling boundaries are fixed at k0, so geometric shrinkage is not imposed by
rescaling the diagnostic domain. Finite-grid enstrophy moments describe this
cylinder only, not a globally identified vortex core or global finite energy.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from meridional_state_cache import StatefulMean, value_and_pressure_modes
from broad_shear_dynamic_control import load_saved_field
from grouped_joined_field import install_in_field


def fixed_cylinder(inner, k0, radial_order=6, axial_order=12, angles=3):
    tau = .5*2.**(-k0)
    ri = inner.from_similarity([inner.p.X_max], [0.], tau)[0, 0]
    zmax = abs(inner.from_similarity([inner.p.X_max], [.5], tau)[0, 2])
    gr, wr = leggauss(radial_order); gz, wz = leggauss(axial_order)
    edges = ri*np.array([0., .1, .2, .4, .8, 1., 2., 4., 8., 16.])
    points, weights = [], []
    for z, zw in zip(zmax*gz, zmax*wz):
        for lo, hi in zip(edges[:-1], edges[1:]):
            for r, rw in zip((lo+hi)/2+(hi-lo)*gr/2, (hi-lo)*wr/2):
                for theta in .31+2*np.pi*np.arange(angles)/angles:
                    points.append([r*np.cos(theta), r*np.sin(theta), z])
                    weights.append(zw*rw*r*2*np.pi/angles)
    return np.asarray(points), np.asarray(weights), dict(
        radius=float(edges[-1]), axial_half_length=float(zmax), reference_k=k0,
        radial_order_per_panel=radial_order, axial_order=axial_order, angles=angles)


def observe(field, points, weights, k):
    tau = .5*2.**(-k); h = 5e-4*np.sqrt(field.nu*tau)
    u = field.fields(points, tau)[0]
    J = np.empty((len(points), 3, 3))
    for i, axis in enumerate(np.eye(3)):
        um2, um, up, up2 = [field.fields(points+j*h*axis, tau)[0]
                            for j in (-2, -1, 1, 2)]
        J[:, :, i] = (um2-8*um+8*up-up2)/(12*h)
    omega = np.column_stack((J[:, 2, 1]-J[:, 1, 2],
                             J[:, 0, 2]-J[:, 2, 0], J[:, 1, 0]-J[:, 0, 1]))
    density = np.sum(omega**2, axis=1)
    enstrophy = weights @ density
    if not np.isfinite(enstrophy) or enstrophy <= 0:
        raise ValueError('Nonpositive or nonfinite sampled enstrophy')
    r = np.linalg.norm(points[:, :2], axis=1)
    swirl = (-points[:, 1]*u[:, 0]+points[:, 0]*u[:, 1])/r
    zmean = weights @ (density*points[:, 2])/enstrophy
    radial_rms = np.sqrt(weights @ (density*r*r)/enstrophy)
    axial_rms = np.sqrt(weights @ (density*(points[:, 2]-zmean)**2)/enstrophy)
    return dict(k=float(k), tau=float(tau), point_count=len(points),
        sampled_vorticity_max=float(np.sqrt(density.max())),
        cylinder_enstrophy=float(enstrophy),
        enstrophy_radial_rms=float(radial_rms),
        enstrophy_axial_rms=float(axial_rms),
        enstrophy_aspect_ratio=float(axial_rms/radial_rms),
        enstrophy_weighted_angular_speed=float(weights @ (density*swirl/r)/enstrophy),
        sampled_swirl_max=float(np.max(abs(swirl))),
        sampled_circulation_max=float(np.max(abs(2*np.pi*r*swirl))),
        cylinder_kinetic_energy=float(.5*(weights @ np.sum(u*u, axis=1))))


def run(input_name='meridional_state_evolution.json', output_name='vortex_state_observables.json'):
    root = Path(__file__).resolve().parent
    evolution = json.loads((root/input_name).read_text())
    dynamic, source = load_saved_field(); install_in_field(dynamic)
    base, values, pressures, _, _ = value_and_pressure_modes(dynamic, source)
    points, weights, domain = fixed_cylinder(base.inner, source['k'])
    rows = []
    # These observables use only velocity. Compact baseline pressure and all
    # additional pressure directions have identically zero velocity, so omit
    # their costly pressure integrals without changing any measured field.
    initial = StatefulMean(dynamic, values, [], np.zeros(len(values)),
                           np.zeros(len(values)), [], source['k'])
    rows.append(dict(path='initial', **observe(initial, points, weights, source['k'])))
    for run in evolution.get('delta_runs', {}).values():
        if run.get('status') != 'completed':
            continue
        for label in ('fixed', 'refreshed'):
            field = StatefulMean(dynamic, values, [], run[label+'_state'],
                run[label+'_slope'], [], run['endpoint_k'])
            rows.append(dict(path=label, delta_k=run['delta_k'],
                             **observe(field, points, weights, run['endpoint_k'])))
    endpoint = evolution.get('endpoint_fit', {})
    if endpoint.get('status') == 'fit_complete':
        field = StatefulMean(dynamic, values, [], endpoint['state'],
                             endpoint['slope'], [], endpoint['k'])
        rows.append(dict(path='constrained', delta_k=evolution['delta_k'],
                         **observe(field, points, weights, endpoint['k'])))
    report = dict(accepted=False, scale_recursion_established=False,
        source_report=input_name, source_status=evolution['status'], domain=domain, rows=rows,
        scope='Fixed physical cylinder, sampled vorticity/enstrophy shape and swirl diagnostics. Not identified core radii, global energy, grid convergence, PDE or scale recursion proof.')
    root.joinpath(output_name).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='meridional_state_evolution.json')
    parser.add_argument('--output', default='vortex_state_observables.json')
    args = parser.parse_args()
    run(args.input, args.output)
