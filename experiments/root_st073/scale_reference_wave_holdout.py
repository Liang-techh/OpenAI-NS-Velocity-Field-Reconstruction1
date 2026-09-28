"""New support-split spatial grid for the refined wave candidate.

Fresh parent field jets plus analytic wave/pressure increments. One sampled
grid is a transfer diagnostic, not a continuum or dynamical certificate.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from grouped_joined_field import install_in_field
from scale_reference_candidate import load_reference
from scale_reference_actual_replay import _spatial_jets, _metric_state, _metric_update, _metric_finish
from scale_reference_refined_wave import _pressure_columns
from supported_fourier_analytic_jets import basis_jets
from scale_transport_generator import solenoidal_scale_generator
from scale_generator_linearization import velocity_correction_response

ROOT = Path(__file__).resolve().parent


def spatial_grid(geometry, order):
    drift = json.loads((ROOT/'localized_drift_1800_fit.json').read_text())['inputs']
    center, widths = np.array(geometry['center']), np.array(geometry['widths'])
    intervals = []
    for axis in range(2):
        low, high = center[axis]-widths[axis], center[axis]+widths[axis]
        edges = [low, high]
        for name in ('first', 'second'):
            c = drift[name+'_patch_center'][axis]
            w = drift[name+'_patch_widths'][axis]
            edges += [max(low, c-w), min(high, c+w)]
        edges = sorted(set(edges))
        nodes, weights = np.polynomial.legendre.leggauss(order)
        intervals.append((np.concatenate([(a+b)/2+(b-a)*nodes/2 for a, b in zip(edges[:-1], edges[1:])]),
                          np.concatenate([(b-a)*weights/2 for a, b in zip(edges[:-1], edges[1:])]), edges))
    r, wr, re = intervals[0]
    z, wz, ze = intervals[1]
    angles = 0.231 + 2*np.pi*np.arange(12)/12
    rr, zz, tt = np.meshgrid(r, z, angles, indexing='ij')
    points = np.column_stack(((rr*np.cos(tt)).ravel(), (rr*np.sin(tt)).ravel(), zz.ravel()))
    weights = (wr[:, None, None]*wz[None, :, None]*rr*2*np.pi/12).ravel()
    return points, weights, dict(order_per_interval=order, radial_edges=re, axial_edges=ze, angular_shift=0.231, angles=12)


def run(order):
    source = ROOT/'scale_reference_refined_wave.json'
    report = json.loads(source.read_text())
    selected = report['selected']
    if selected is None:
        raise ValueError('No wave candidate')
    geometry = json.loads((ROOT/'full_wave_frozen_cache.json').read_text())['inputs']['wave']
    reference = load_reference(step=ROOT/'scale_reference_trust_nonlinear_fit.json')
    install_in_field(reference.base_field)
    packed = np.asarray(geometry['coefficients_original'])
    coefficients = packed[:, 0] + 1j*packed[:, 1]
    amplitude = selected['amplitude_shift']
    pressure = np.asarray(selected['pressure_coefficients'])
    points, weights, definition = spatial_grid(geometry, order)
    states = {key: _metric_state() for key in ('parent', 'candidate')}
    for begin in range(0, len(points), 768):
        p, w = points[begin:begin+768], weights[begin:begin+768]
        u, j, gp, lap = _spatial_jets(reference, p, reference.tau0, 1e-6)
        theta = np.arctan2(p[:, 1], p[:, 0])
        swirl = (-np.sin(theta)*u[:, 0]+np.cos(theta)*u[:, 1]).reshape(-1, 12).mean(axis=1).repeat(12)
        residual = solenoidal_scale_generator(p, u, j, reference.tau0, reference.h, swirl)
        residual += np.einsum('nij,nj->ni', j, u)-reference.nu*lap+gp
        v, dj, visc, _, _ = basis_jets(p, geometry['center'], geometry['widths'], 1, geometry['degree'], geometry['carrier'], reference.nu)
        v = amplitude*np.einsum('niq,q->ni', v, coefficients).real
        dj = amplitude*np.einsum('nijq,q->nij', dj, coefficients).real
        lv = -amplitude*np.einsum('niq,q->ni', visc, coefficients).real/reference.nu
        pressure_design, _ = _pressure_columns(p, np.asarray(geometry['carrier']))
        linear, quadratic = velocity_correction_response(p, u, j, v, dj, lv,
            reference.tau0, reference.h, np.zeros(len(p)), reference.nu,
            (pressure_design@pressure).reshape(-1, 3))
        _metric_update(states['parent'], residual, w, p)
        _metric_update(states['candidate'], residual+linear+quadratic, w, p)
        print(f'{min(begin+768, len(points))}/{len(points)}', flush=True)
    metrics = {key: _metric_finish(value, len(points), float(weights.sum())) for key, value in states.items()}
    result = dict(status='completed', accepted=False, pde_validated=False, scale_recursion_established=False,
        candidate_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), grid=definition,
        metrics=metrics, L2_ratio=metrics['candidate']['volume_L2']/metrics['parent']['volume_L2'],
        max_ratio=metrics['candidate']['max_norm']/metrics['parent']['max_norm'],
        scope='New support-split grid over the same full wave cylinder. Parent spatial FD and analytic increments; zero force, scale-generator time derivative. Spatial convergence and actual candidate FD agreement remain separate requirements.')
    (ROOT/f'scale_reference_wave_holdout_o{order}.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--order', type=int, default=4)
    run(parser.parse_args().order)
