"""Resolved angular momentum-defect budget from the actual reference jet cache."""
import hashlib
import json
from pathlib import Path

import numpy as np
from wave_residual_harmonics import budget

ROOT = Path(__file__).resolve().parent


def run():
    source = ROOT / 'scale_reference_trust_nonlinear_actual.json'
    replay = json.loads(source.read_text())
    cache = replay['chunk_caches']['corrected_reference']
    directory = ROOT / '.scale_replay_cache' / cache['fingerprint']
    keys = ('points', 'weights', 'velocity', 'jacobian', 'laplacian',
            'pressure_gradient', 'generator', 'residual')
    pieces = {key: [] for key in keys}
    manifest = []
    count = 0
    for path in sorted(directory.glob('*.npz')):
        if int(path.stem) != count:
            raise ValueError('Noncontiguous chunks')
        with np.load(path, allow_pickle=False) as data:
            for key in keys:
                pieces[key].append(data[key])
            count += len(data['points'])
        manifest.append(dict(name=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    data = {key: np.concatenate(value) for key, value in pieces.items()}
    p, w, u, j, residual = (data[key] for key in ('points', 'weights', 'velocity', 'jacobian', 'residual'))
    angular = budget(p, w, residual, 12)  # Also validates full angular rings.
    expected = replay['cases']['corrected_reference']['metric']
    if not np.isclose(angular['momentum_volume_L2'], expected['volume_L2'], rtol=1e-12):
        raise ValueError('Cached residual does not reproduce frozen report')
    theta = np.arctan2(p[:, 1], p[:, 0])
    er = np.column_stack((np.cos(theta), np.sin(theta), np.zeros(count)))
    et = np.column_stack((-np.sin(theta), np.cos(theta), np.zeros(count)))
    def mean(values):
        return values.reshape(-1, 12).mean(axis=1)
    def azimuthal(values):
        return mean(np.sum(et * values, axis=1))
    rw = w.reshape(-1, 12).sum(axis=1)
    target = azimuthal(residual)
    square = float(np.sum(rw * target**2))
    def metric(values):
        return dict(volume_L2=float(np.sqrt(np.sum(rw * values**2))),
                    signed_projection=float(np.sum(rw * values * target) / square))
    nu = replay['inputs']['viscosity']
    convection = np.einsum('nij,nj->ni', j, u)
    terms = dict(generator=data['generator'], convection=convection,
                 viscosity=-nu*data['laplacian'], pressure=data['pressure_gradient'])
    ur, uth, uz = mean(np.sum(er*u, axis=1)), azimuthal(u), mean(u[:, 2])
    dr = mean(np.einsum('ni,nij,nj->n', et, j, er))
    dz = mean(np.sum(et*j[:, :, 2], axis=1))
    radii = np.hypot(p[:, 0], p[:, 1]).reshape(-1, 12)[:, 0]
    mean_convection = ur*dr + uz*dz + ur*uth/radii
    report = dict(status='completed', accepted=False, pde_validated=False,
                  scale_recursion_established=False,
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  cache_fingerprint=cache['fingerprint'], chunk_manifest=manifest,
                  angular=angular,
                  mean_azimuthal=dict(volume_L2=float(np.sqrt(square)),
                      squared_L2_fraction=square/angular['momentum_volume_L2']**2),
                  mean_azimuthal_terms={name: metric(azimuthal(value)) for name, value in terms.items()},
                  reynolds_split=dict(mean_flow=metric(mean_convection),
                      fluctuations=metric(azimuthal(convection)-mean_convection)),
                  scope='Sampled corrected-field budget. Angular mean azimuthal pressure gradient vanishes analytically for periodic single-valued pressure; numerical norms are not continuum lower bounds. Term norms are not additive.')
    output = ROOT / 'scale_reference_refined_budget.json'
    output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: report[key] for key in ('mean_azimuthal', 'mean_azimuthal_terms', 'reynolds_split')}))
    return report


if __name__ == '__main__':
    run()
