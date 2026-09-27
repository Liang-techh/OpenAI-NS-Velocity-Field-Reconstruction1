"""Whole-support energy quadrature for the localized candidate, not PDE acceptance."""
import hashlib
import json
import os
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from global_axial_extension import build_candidate

ROOT = Path(__file__).resolve().parent


def support_bounds(localized, geometry, tau, extra_boxes=()):
    """Bound the moving mean, fixed wave, and any added (r0,r1,z0,z1) boxes."""
    if not localized.tau_min <= tau <= localized.tau_max:
        raise ValueError('Time outside registered field domain')
    eta = localized.eta_outer
    # Obtain physical z from the same coordinate map as the actual inner field.
    point = localized.inner.from_similarity([localized.join_X], [eta], tau)[0]
    zmean = abs(float(point[2]))
    rmean = localized.radial_radii(tau)[1]
    center, width = np.asarray(geometry['center']), np.asarray(geometry['widths'])
    extra_boxes = np.asarray(extra_boxes, dtype=float).reshape(-1,4)
    if len(extra_boxes) and (not np.isfinite(extra_boxes).all() or
            np.any(extra_boxes[:,0]<0) or np.any(extra_boxes[:,1]<=extra_boxes[:,0]) or
            np.any(extra_boxes[:,3]<=extra_boxes[:,2])):
        raise ValueError('Invalid additional support boxes')
    return dict(radius=max([rmean, float(center[0]+width[0])] + extra_boxes[:,1].tolist()),
                z_lower=min([-zmean, float(center[1]-width[1])] + extra_boxes[:,2].tolist()),
                z_upper=max([zmean, float(center[1]+width[1])] + extra_boxes[:,3].tolist()),
                mean_radius=rmean, mean_z_halfwidth=zmean)


def panel_nodes(edges, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x, w = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        x.extend((a+b)/2 + (b-a)/2*nodes)
        w.extend((b-a)/2*weights)
    return np.asarray(x), np.asarray(w)


def run():
    full, localized, _, candidate, snapshot, _ = build_candidate()
    tau = snapshot['inputs']['mean']['tau']
    geometry = snapshot['inputs']['wave']
    bounds = support_bounds(localized, geometry, tau)
    center, width = np.asarray(geometry['center']), np.asarray(geometry['widths'])
    rflat = localized.radial_radii(tau)[0]
    redges = sorted(set([0., center[0]-width[0], center[0],
                        center[0]+width[0], rflat, bounds['radius']]))
    zedges = sorted(set([bounds['z_lower'], center[1]-width[1], 0.,
                        center[1]+width[1], bounds['z_upper']]))
    report = dict(status='running', accepted=False, pde_validated=False,
                  scale_recursion_established=False, tau=tau, bounds=bounds,
                  scope='Whole bounding-cylinder kinetic energy quadrature at one time. No certified integration error, uniform terminal energy bound, or NS acceptance.',
                  source_hashes={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                 for name in ('global_axial_extension.py', 'enriched_mean_endpoint_tangent.json')},
                  rows=[])
    output = ROOT/'global_support_energy.json'
    for order in (6, 10):
        radius, rw = panel_nodes(redges, order)
        z, zw = panel_nodes(zedges, order)
        angles = np.arange(8)*2*np.pi/8
        rr, zz, aa = np.meshgrid(radius, z, angles, indexing='ij')
        points = np.column_stack((rr.ravel()*np.cos(aa.ravel()),
                                  rr.ravel()*np.sin(aa.ravel()), zz.ravel()))
        weights = (rw[:,None,None]*zw[None,:,None]*np.ones((1,1,8))*rr*2*np.pi/8).ravel()
        velocity = full.fields(points,tau)[0]
        if not np.isfinite(velocity).all():
            raise ValueError('Nonfinite velocity in support integration')
        energy = .5*float(weights @ np.sum(velocity*velocity,axis=1))
        report['rows'].append(dict(order=order, points=len(points), kinetic_energy=energy,
                                   sampled_speed_max=float(np.max(np.linalg.norm(velocity,axis=1)))))
        output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report['rows'][-1]),flush=True)
    report['relative_energy_difference'] = abs(report['rows'][1]['kinetic_energy']-report['rows'][0]['kinetic_energy'])/report['rows'][1]['kinetic_energy']
    report['status']='completed'
    output.write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    run()
