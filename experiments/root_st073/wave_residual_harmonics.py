"""Physical-volume Fourier budget of a sampled Cartesian momentum residual.

Rotate vectors into cylindrical components before decomposing angular modes.
This is postprocessing, not a residual or continuum accuracy certificate.
Samples must be contiguous, equally spaced angular rings with equal weights.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def budget(points, weights, residual, angles):
    points = np.asarray(points, float)
    weights = np.asarray(weights, float)
    residual = np.asarray(residual, float)
    if points.shape != residual.shape or points.ndim != 2 or points.shape[1] != 3:
        raise ValueError('points and residual must have shape (N,3)')
    if len(points) % angles or weights.shape != (len(points),):
        raise ValueError('Incomplete angular rings or invalid weights')
    if angles < 3 or not np.all(np.isfinite(residual)) or np.any(weights <= 0):
        raise ValueError('Require finite residual, positive weights and >=3 angles')
    rings = points.reshape(-1, angles, 3)
    wr = weights.reshape(-1, angles)
    radius = np.hypot(rings[..., 0], rings[..., 1])
    theta = np.arctan2(rings[..., 1], rings[..., 0])
    if np.any(radius <= 0):
        raise ValueError('Cylindrical decomposition requires positive radius')
    if not np.allclose(radius, radius[:, :1], rtol=1e-10, atol=1e-14):
        raise ValueError('Ring radii differ')
    if not np.allclose(rings[..., 2], rings[:, :1, 2], rtol=1e-10, atol=1e-14):
        raise ValueError('Ring z coordinates differ')
    phase = np.exp(1j * (theta - theta[:, :1]))
    if not np.allclose(phase, np.exp(2j*np.pi*np.arange(angles)/angles), atol=1e-10):
        raise ValueError('Angles must increase uniformly within each ring')
    if not np.allclose(wr, wr[:, :1], rtol=1e-10, atol=0):
        raise ValueError('Angular weights must be equal within each ring')
    rr = residual.reshape(-1, angles, 3)
    c, s = np.cos(theta), np.sin(theta)
    cylindrical = np.stack((c*rr[...,0]+s*rr[...,1],
                            -s*rr[...,0]+c*rr[...,1], rr[...,2]), axis=-1)
    coefficients = np.fft.rfft(cylindrical, axis=1)/angles
    ring_weights = wr.sum(axis=1)
    total = float(np.sum(weights[:,None]*residual**2))
    rows = []
    for mode in range(coefficients.shape[1]):
        multiplicity = 1 if mode == 0 or (angles % 2 == 0 and mode == angles//2) else 2
        components = multiplicity*np.sum(ring_weights[:,None]*abs(coefficients[:,mode])**2, axis=0)
        square = float(components.sum())
        rows.append(dict(mode=mode, volume_L2=float(np.sqrt(square)),
                         squared_L2_fraction=square/total if total else 0.0,
                         cylindrical_component_squared_L2=components.tolist()))
    parsed = sum(row['volume_L2']**2 for row in rows)
    return dict(point_count=len(points), angles=angles, physical_volume=float(weights.sum()),
                momentum_max=float(np.linalg.norm(residual,axis=1).max()),
                momentum_volume_L2=float(np.sqrt(total)), modes=rows,
                parseval_relative_error=abs(parsed-total)/max(total,1e-300),
                scope='Sampled cylindrical Fourier budget; angular aliasing and spatial accuracy unverified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cache', type=Path)
    parser.add_argument('--prefix', default='hold')
    parser.add_argument('--angles', type=int, default=12)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with np.load(args.cache, allow_pickle=False) as data:
        result = budget(data[args.prefix+'_points'], data[args.prefix+'_weights'],
                        data[args.prefix+'_residual'], args.angles)
    result['source_cache'] = str(args.cache)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
