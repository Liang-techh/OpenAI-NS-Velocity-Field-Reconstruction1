"""Replay harmonic and SVD-truncation diagnostics from frozen physical jets.

No additional physical jets are evaluated. Predicted tangent residuals are
not substitutes for independent actual-field finite-difference replay.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from supported_fourier_basis import basis_data
from wave_residual_harmonics import budget

ROOT = Path(__file__).resolve().parent


def _basis_columns(points, center, widths, carriers, degree):
    blocks = []
    for mode in (0,1,2):
        velocity,_,gradient = basis_data(points,center,widths,mode,degree,carriers[mode])
        for tensor in (velocity,gradient):
            for j in range(tensor.shape[-1]):
                blocks.append(tensor[:,:,j].real.reshape(-1))
                if mode:
                    blocks.append(-tensor[:,:,j].imag.reshape(-1))
    return np.stack(blocks,axis=1), None


def run():
    raw = (ROOT / 'full_wave_frozen_cache.json').read_bytes()
    source = json.loads(raw)
    cache = source['frozen_cache']
    wave = source['inputs']['wave']
    carriers = {0: np.zeros(2), 1: np.array(wave['carrier']),
                2: 2*np.array(wave['carrier'])}
    designs = {}
    residuals = {}
    weights = {}
    for label in ('fit', 'holdout'):
        designs[label], _ = _basis_columns(np.array(cache[label+'_points']),
            wave['center'], wave['widths'], carriers, wave['degree'])
        residuals[label] = np.array(cache[label+'_residual'])
        weights[label] = np.array(cache[label+'_weights'])
    weighted = designs['fit']*np.repeat(np.sqrt(weights['fit']),3)[:,None]
    scales = np.maximum(np.linalg.norm(weighted,axis=0),1e-30)
    u,s,v = np.linalg.svd(weighted/scales,full_matrices=False)
    rhs = (-residuals['fit']*np.sqrt(weights['fit'][:,None])).ravel()
    result = dict(accepted=False,pde_validated=False,scale_recursion_established=False,
                  source_sha256=hashlib.sha256(raw).hexdigest(),inputs=source['inputs'],
                  singular_values=s.tolist(),column_count=weighted.shape[1],
                  scope='Frozen physical FD residual and predicted linear tangent correction; actual-field replay is separate.',
                  truncation_scope='Same holdout used for screening; not independent validation of selected truncation.')
    rows = []
    for tol in (1e-10,1e-8,1e-6,1e-4,1e-3,1e-2,.05,.1):
        keep = s > tol*s[0]
        x = (v[keep].T@((u[:,keep].T@rhs)/s[keep]))/scales
        row = dict(rcond=tol,rank=int(keep.sum()))
        for label, angles in (('fit',8),('holdout',12)):
            corrected = residuals[label]+(designs[label]@x).reshape(-1,3)
            metrics = budget(cache[label+'_points'],weights[label],corrected,angles)
            row[label+'_L2'] = metrics['momentum_volume_L2']
            row[label+'_max'] = metrics['momentum_max']
            if tol == 1e-10:
                result[label] = dict(frozen=budget(cache[label+'_points'],weights[label],
                    residuals[label],angles),predicted_full180=metrics)
        if tol == 1e-10:
            result.update(rank=int(keep.sum()),coefficients=x.tolist())
        rows.append(row)
    result['truncation_diagnostic'] = rows
    (ROOT/'full_wave_harmonic_budget.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(rows,indent=2))


if __name__ == '__main__':
    run()
