"""Independent actual-field replay of the momentum-aware wave co-design."""
import hashlib
import argparse
import json
import time
from pathlib import Path

import numpy as np

from affine_momentum import jets,momentum
from broad_meridional_constrained import load_saved_field
from full_wave_tangent import LocalPotentialField,_unpack_full
from grouped_joined_field import install_in_field
from wave_residual_harmonics import budget

ROOT = Path(__file__).resolve().parent


def run(source_path=None,output_path=None):
    start = time.perf_counter()
    source_path = Path(source_path) if source_path else ROOT/'wave_dynamics_codesign.json'
    source_raw = source_path.read_bytes()
    source = json.loads(source_raw)
    if source['status'] != 'completed':
        raise ValueError('Finish and freeze optimization before independent replay')
    dense = json.loads((ROOT/'full_wave_dense_tangent.json').read_text())
    if dense['status'] != 'completed':
        raise ValueError('Independent dense grid must be complete')
    snapshot = json.loads((ROOT/'full_wave_frozen_cache.json').read_text())
    wave = snapshot['inputs']['wave']
    mean,mean_report = load_saved_field()
    install_in_field(mean)
    if mean_report['coefficients'] != snapshot['inputs']['mean']['coefficients']:
        raise ValueError('Mean coefficients changed since training cache')
    raw = np.array(source['selected']['coefficients_original'])
    coefficient = raw[:,0]+1j*raw[:,1]
    q = (wave['degree']+1)**2
    initial = (np.zeros(3*q,complex),coefficient,np.zeros(3*q,complex))
    derivative,pressure = _unpack_full(source['selected']['tangent_coefficients'],q)
    carrier = np.array(wave['carrier'])
    tau = snapshot['inputs']['mean']['tau']
    field = LocalPotentialField(mean,wave['center'],wave['widths'],wave['degree'],
        {0:np.zeros(2),1:carrier,2:2*carrier},initial,derivative,pressure,tau)
    data = dense['new_frozen_cache']
    points,weights = np.array(data['points']),np.array(data['weights'])
    report = dict(status='running',accepted=False,pde_validated=False,
        scale_recursion_established=False,constraints_maintained=False,
        source_sha256=hashlib.sha256(source_raw).hexdigest(),
        source=source_path.name,inputs=snapshot['inputs'],
        selected_coefficients_original=source['selected']['coefficients_original'],
        tangent_coefficients=source['selected']['tangent_coefficients'],
        independent_geometry=dense['new_grid_definition'],
        previous_candidate_frozen=dense['new_harmonic_budget_frozen'],
        previous_candidate_corrected=dense['new_grid_actual_selected']['metric'],
        scope='Actual Cartesian full momentum at one time on a spatial grid unused by this shape optimization. No continuous trajectory, moment/cone compatibility, or recursion acceptance.')
    output = Path(output_path) if output_path else ROOT/'wave_dynamics_replay.json'
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Independent complete-field replay started',flush=True)
    residual = momentum(jets(field,points,tau,snapshot['timesteps']['hspace'],
                              snapshot['timesteps']['htime']))
    report['momentum'] = budget(points,weights,residual,12)
    report.update(status='completed',elapsed_seconds=time.perf_counter()-start)
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report['momentum'],indent=2),flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source',type=Path)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    run(args.source,args.output)
