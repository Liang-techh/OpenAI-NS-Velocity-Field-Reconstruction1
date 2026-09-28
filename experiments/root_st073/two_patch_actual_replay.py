"""Actual independent momentum replay for a frozen drift-capped two-patch fit."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from affine_momentum import jets, momentum
from enriched_shape_replay import build_field
from full_wave_tangent import LocalPotentialField, _unpack_full
from localized_actual_replay import metric

ROOT=Path(__file__).resolve().parent


def load():
    path=ROOT/'localized_two_patch_constrained.json'
    raw=path.read_bytes();report=json.loads(raw)
    if report['status']!='completed' or not report['assembled_feasible'] or not report['inputs']['drift_cap_enabled']:
        raise ValueError('Need a frozen feasible drift-capped two-patch candidate')
    parent_path=ROOT/report['sources']['candidate']['path']
    if hashlib.sha256(parent_path.read_bytes()).hexdigest()!=report['sources']['candidate']['sha256']:
        raise ValueError('Parent mismatch')
    field,snapshot=build_field(json.loads(parent_path.read_text()))
    tau0=snapshot['inputs']['mean']['tau']
    carrier=np.asarray(snapshot['inputs']['wave']['carrier'])
    for number,word in ((1,'first'),(2,'second')):
        derivative,pressure=_unpack_full(np.asarray(report['selected'][f'patch{number}_coefficients']),9)
        field=LocalPotentialField(field,report['inputs'][f'{word}_patch_center'],
                                  report['inputs'][f'{word}_patch_widths'],2,
                                  {0:np.zeros(2),1:carrier,2:2*carrier},
                                  tuple(np.zeros(27,complex) for _ in range(3)),
                                  derivative,pressure,tau0)
    return field,snapshot,report,hashlib.sha256(raw).hexdigest()


def run():
    started=time.perf_counter()
    field,snapshot,candidate,source_hash=load()
    baseline_path=ROOT/'localized_actual_replay.json'
    baseline=json.loads(baseline_path.read_text())
    grid_path=ROOT/'full_wave_dense_tangent.json'
    grid_hash=hashlib.sha256(grid_path.read_bytes()).hexdigest()
    if (baseline['status']!='completed' or baseline['grid_sha256']!=grid_hash
        or baseline['parent_sha256']!=candidate['sources']['candidate']['sha256']
        or baseline['source_sha256']!=candidate['sources']['first_fit']['sha256']
        or baseline['timesteps']!=snapshot['timesteps']):
        raise ValueError('Frozen actual baseline is incompatible')
    grid=json.loads(grid_path.read_text())['new_frozen_cache']
    points,weights=np.asarray(grid['points']),np.asarray(grid['weights'])
    output=ROOT/'two_patch_actual_replay.json'
    report=dict(status='running',accepted=False,pde_validated=False,scale_recursion_established=False,
                source_sha256=source_hash,baseline_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
                grid_sha256=grid_hash,point_count=len(points),rows=[],
                scope='Actual five-point Cartesian finite differences at two times on independent grid. Baseline metrics reused only after parent/grid/first-fit/timestep hash checks; no full-support or interval certificate.')
    def save():output.write_text(json.dumps(report,indent=2)+'\n')
    save()
    for previous in baseline['rows']:
        dk=previous['delta_k'];tau=snapshot['inputs']['mean']['tau']*2**(-dk)
        if tau!=previous['tau']:raise ValueError('Baseline evaluation time differs')
        print(json.dumps(dict(stage='running',delta_k=dk)),flush=True)
        residual=momentum(jets(field,points,tau,snapshot['timesteps']['hspace'],snapshot['timesteps']['htime']))
        result=metric(residual,weights)
        row=dict(delta_k=dk,tau=tau,parent=previous['parent'],first_patch=previous['corrected'],two_patch=result,
                 both_metrics_improve_vs_first=all(result[k]<previous['corrected'][k] for k in result))
        report['rows'].append(row);save();print(json.dumps(row),flush=True)
    report.update(status='completed',elapsed_seconds=time.perf_counter()-started);save()


if __name__=='__main__':run()
