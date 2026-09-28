"""Add the exact quadratic convection remainder to cached endpoint momentum.

This uses the same fitted grid and cached finite-difference parent derivatives.
It is not an independent-grid or continuum residual certificate.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from pathlib import Path
import numpy as np
from endpoint_acceleration_projection import _unpack_acceleration_control
from supported_fourier_analytic_jets import basis_jets

ROOT = Path(__file__).resolve().parent


def metric(residual, weights):
    squared = np.sum(residual * residual, axis=1)
    return {'maximum': float(np.sqrt(squared.max())),
            'volume_L2': float(np.sqrt(weights @ squared))}


def run():
    source = ROOT / 'shape_constrained_acceleration.json'
    candidate = json.loads(source.read_text())
    if candidate['status'] != 'completed':
        raise ValueError('Wait for a frozen completed acceleration solve')
    control = np.asarray(candidate['selected']['acceleration_coefficients'])
    geometry_path = ROOT / 'full_wave_frozen_cache.json'
    geometry = json.loads(geometry_path.read_text())['inputs']['wave']
    projection_path = ROOT / 'endpoint_acceleration_projection.json'
    projection = json.loads(projection_path.read_text())
    dt = projection['inputs']['physical_time_increment']
    nu = projection['inputs']['viscosity']
    cache_path = ROOT / 'acceleration_momentum_cache.npz'
    with np.load(cache_path, allow_pickle=False) as cache:
        if str(cache['design_convention']) != 'real-minus-imag':
            raise ValueError('Stale complex column convention')
        if str(cache['projection_sha256']) != hashlib.sha256(projection_path.read_bytes()).hexdigest():
            raise ValueError('Projection source changed after cache construction')
        points, weights = cache['points'], cache['weights']
        base = cache['base_residual']
        linear = base + (cache['design'] @ control).reshape(-1, 3)
    unpacked = _unpack_acceleration_control(control)
    du = np.zeros_like(points)
    dJ = np.zeros((len(points), 3, 3))
    for mode, (acceleration, _pressure) in unpacked.items():
        if not np.any(acceleration):
            continue
        V, J, _, _, _ = basis_jets(points, geometry['center'], geometry['widths'],
                                   mode, 2, mode * np.asarray(geometry['carrier']), nu)
        du += 0.5 * dt**2 * np.real(np.einsum('ncq,q->nc', V, acceleration))
        dJ += 0.5 * dt**2 * np.real(np.einsum('ncdq,q->ncd', J, acceleration))
    remainder = np.einsum('nij,nj->ni', dJ, du)
    nonlinear = linear + remainder
    cap = candidate['selected']['momentum_peak_cap']
    report = dict(status='completed', accepted=False, pde_validated=False,
                  scale_recursion_established=False,
                  scope='Same fitted grid; cached FD parent plus analytic correction and exact quadratic convection. No independent spatial/time replay.',
                  sources={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in (source, geometry_path, projection_path, cache_path)},
                  point_count=len(points), base=metric(base, weights),
                  linearized=metric(linear, weights), nonlinear=metric(nonlinear, weights),
                  quadratic_remainder=metric(remainder, weights), peak_cap=cap,
                  nonlinear_peak_cap_pass=bool(np.linalg.norm(nonlinear, axis=1).max() <= cap * (1 + 1e-10)),
                  velocity_correction_max=float(np.linalg.norm(du, axis=1).max()),
                  shape_margins=candidate['selected']['shape_margin'])
    (ROOT / 'acceleration_nonlinear_replay.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
