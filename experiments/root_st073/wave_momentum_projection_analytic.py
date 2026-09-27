"""Replace wave finite differences by analytic derivatives, reusing mean jets.

The mean remains finite-differenced. This is not a complete analytic NS
residual, and agreement with the older cache is not a PDE tolerance test.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from supported_fourier_analytic_jets import basis_jets
from wave_momentum_projection import _assemble_wave_columns,residual_and_jacobian,_metric

ROOT = Path(__file__).resolve().parent


def run():
    source = ROOT/'wave_momentum_projection.npz'
    metadata_raw = (ROOT/'wave_momentum_projection.json').read_bytes()
    metadata = json.loads(metadata_raw)
    if metadata['status'] != 'completed':
        raise ValueError('Finish the source projection cache first')
    with np.load(source,allow_pickle=False) as data:
        cache = {key:data[key] for key in data.files}
    wave = json.loads((ROOT/'wave_stress_growth_codesign.json').read_text())
    c = np.array(wave['selected']['coefficients_original'])
    x = np.r_[c[:,0],c[:,1]]
    previous,_ = residual_and_jacobian(cache,x)
    A,B,L,layout = _assemble_wave_columns(cache['points'],cache['mean_velocity'],
        cache['mean_gradient'],wave['center'],wave['widths'],wave['carrier'],
        wave['degree'],wave['viscosity'],metadata['timesteps']['hspace'],
        jet_function=basis_jets)
    cache.update(A=A,B=B,L=L)
    current,_ = residual_and_jacobian(cache,x)
    output_cache = ROOT/'wave_momentum_projection_analytic.npz'
    np.savez_compressed(output_cache,**cache)
    report = dict(status='completed',accepted=False,pde_validated=False,
        scale_recursion_established=False,
        source_metadata_sha256=hashlib.sha256(metadata_raw).hexdigest(),
        source_cache_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        inputs=metadata['inputs'],backend='analytic compact wave derivatives; cached FD mean',
        cache_path=output_cache.name,layout=layout,
        frozen_previous=_metric(previous,cache['weights']),
        frozen_analytic_wave=_metric(current,cache['weights']),
        difference=_metric(current-previous,cache['weights']),
        scope='Same instantaneous original wave, with analytic wave derivatives and unchanged mean finite differences. Tangent velocity/pressure columns are unchanged. No full analytic or PDE certificate.')
    (ROOT/'wave_momentum_projection_analytic.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report['difference'],indent=2))
    return report


if __name__ == '__main__':
    run()
