"""Short-interval anisotropic profile drift, not a recursion certificate."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import json
import hashlib
from pathlib import Path
import numpy as np
from global_localized_candidate import load

ROOT = Path(__file__).resolve().parent


def cylindrical(velocity, points):
    theta = np.arctan2(points[:, 1], points[:, 0])
    c, s = np.cos(theta), np.sin(theta)
    return np.column_stack((c*velocity[:, 0]+s*velocity[:, 1],
                            -s*velocity[:, 0]+c*velocity[:, 1], velocity[:, 2]))


def run():
    full, parent, localized, snapshot, _, _, hashes = load()
    cache_path = ROOT/'enriched_endpoint_shape_cache.npz'
    with np.load(cache_path, allow_pickle=False) as cache:
        points, weights = cache['points'], cache['weights']
    tau0 = snapshot['inputs']['mean']['tau']
    h = float(localized.inner.h)
    report = dict(status='running', accepted=False, scale_recursion_established=False,
                  pde_validated=False, source_hashes=hashes,
                  grid_sha256=hashlib.sha256(cache_path.read_bytes()).hexdigest(),
                  point_count=len(points), h=h,
                  convention='Fixed similarity labels: xy*=s^.5,z*=s^(.5-h),tau*=s; pull back cylindrical velocity with s^.5,s^(.5+h),s^(.5+h). Reference-volume weights.',
                  caveat='Measures local profile stationarity, not general recursion. Dynamic amplitudes and log-periodic profiles may have nonzero drift. Sampled tiny time intervals only.',
                  rows=[])
    for name, field in (('parent', parent), ('localized', full)):
        reference = cylindrical(field.fields(points, tau0)[0], points)
        norm = np.sqrt(weights @ np.sum(reference**2, axis=1))
        component_norm = np.sqrt(weights @ (reference**2))
        for dk in (1e-6, 2e-6, 1e-5):
            s = 2**(-dk)
            mapped = points * np.array([s**.5, s**.5, s**(.5-h)])
            velocity = cylindrical(field.fields(mapped, tau0*s)[0], mapped)
            pulled = velocity * np.array([s**.5, s**(.5+h), s**(.5+h)])
            error = pulled-reference
            relative = float(np.sqrt(weights @ np.sum(error**2, axis=1))/norm)
            components = np.sqrt(weights @ (error**2))/np.maximum(component_norm, 1e-300)
            row = dict(field=name, delta_k=dk, scale=s, relative_profile_L2=relative,
                       relative_L2_per_abs_log_scale=relative/abs(np.log(s)),
                       component_relative_L2=components.tolist())
            report['rows'].append(row)
            print(json.dumps(row), flush=True)
    report['status']='completed'
    (ROOT/'localized_similarity_drift.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
